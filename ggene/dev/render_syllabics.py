#!/usr/bin/env python3
"""
Render Canadian Aboriginal Syllabics characters as PNG images for visual analysis.
"""

from PIL import Image, ImageDraw, ImageFont
import os
from pathlib import Path

# Characters to organize (from vocab.py)
CHAR_HEX = [i for i in range(0x1400, 0x1680)] + [i for i in range(0x18B0, 0x18F6)]
CHARACTERS = "".join([chr(i) for i in CHAR_HEX])
_CHARACTERS = """
ᐁ ᐂ ᐃ ᐄ ᐅ ᐆ ᐇ ᐈ ᐉ ᐊ ᐋ ᐌ ᐍ ᐎ ᐏ ᐐ ᐑ ᐒ ᐓ ᐔ ᐕ ᐖ ᐗ ᐘ ᐙ ᐚ ᐛ ᐜ ᐝ ᐞ ᐟ
ᐠ ᐡ ᐢ ᐣ ᐤ ᐥ ᐦ ᐧ ᐨ ᐩ ᐪ ᐫ ᐬ ᐭ ᐮ ᐯ ᐰ ᐱ ᐲ ᐳ ᐴ ᐵ ᐶ ᐷ ᐸ ᐹ ᐺ ᐻ ᐼ ᐽ ᐾ ᐿ
ᑀ ᑁ ᑂ ᑃ ᑄ ᑅ ᑆ ᑇ ᑈ ᑉ ᑊ ᑋ ᑌ ᑍ ᑎ ᑏ ᑐ ᑑ ᑒ ᑓ ᑔ ᑕ ᑖ ᑗ ᑘ ᑙ ᑚ ᑛ ᑜ ᑝ ᑞ ᑟ
ᑠ ᑡ ᑢ ᑣ ᑤ ᑥ ᑦ ᑧ ᑨ ᑩ ᑪ ᑫ ᑬ ᑭ ᑮ ᑯ ᑰ ᑱ ᑲ ᑳ ᑴ ᑵ ᑶ ᑷ ᑸ ᑹ ᑺ ᑻ ᑼ ᑽ ᑾ ᑿ
ᒀ ᒁ ᒂ ᒃ ᒄ ᒅ ᒆ ᒇ ᒈ ᒉ ᒊ ᒋ ᒌ ᒍ ᒎ ᒏ ᒐ ᒑ ᒒ ᒓ ᒔ ᒕ ᒖ ᒗ ᒘ ᒙ ᒚ ᒛ ᒜ ᒝ ᒞ ᒟ
ᒠ ᒡ ᒢ ᒣ ᒤ ᒥ ᒦ ᒧ ᒨ ᒩ ᒪ ᒫ ᒬ ᒭ ᒮ ᒯ ᒰ ᒱ ᒲ ᒳ ᒴ ᒵ ᒶ ᒷ ᒸ ᒹ ᒺ ᒻ ᒼ ᒽ ᒾ ᒿ
ᓀ ᓁ ᓂ ᓃ ᓄ ᓅ ᓆ ᓇ ᓈ ᓉ ᓊ ᓋ ᓌ ᓍ ᓎ ᓏ ᓐ ᓑ ᓒ ᓓ ᓔ ᓕ ᓖ ᓗ ᓘ ᓙ ᓚ ᓛ ᓜ ᓝ ᓞ ᓟ
ᓠ ᓡ ᓢ ᓣ ᓤ ᓥ ᓦ ᓧ ᓨ ᓩ ᓪ ᓫ ᓬ ᓭ ᓮ ᓯ ᓰ ᓱ ᓲ ᓳ ᓴ ᓵ ᓶ ᓷ ᓸ ᓹ ᓺ ᓻ ᓼ ᓽ ᓾ ᓿ
ᔀ ᔁ ᔂ ᔃ ᔄ ᔅ ᔆ ᔇ ᔈ ᔉ ᔊ ᔋ ᔌ ᔍ ᔎ ᔏ ᔐ ᔑ ᔒ ᔓ ᔔ ᔕ ᔖ ᔗ ᔘ ᔙ ᔚ ᔛ ᔜ ᔝ ᔞ ᔟ
ᔠ ᔡ ᔢ ᔣ ᔤ ᔥ ᔦ ᔧ ᔨ ᔩ ᔪ ᔫ ᔬ ᔭ ᔮ ᔯ ᔰ ᔱ ᔲ ᔳ ᔴ ᔵ ᔶ ᔷ ᔸ ᔹ ᔺ ᔻ ᔼ ᔽ ᔾ ᔿ
ᕀ ᕁ ᕂ ᕃ ᕄ ᕅ ᕆ ᕇ ᕈ ᕉ ᕊ ᕋ ᕌ ᕍ ᕎ ᕏ ᕐ ᕑ ᕒ ᕓ ᕔ ᕕ ᕖ ᕗ ᕘ ᕙ ᕚ ᕛ ᕜ ᕝ ᕞ ᕟ
ᕠ ᕡ ᕢ ᕣ ᕤ ᕥ ᕦ ᕧ ᕨ ᕩ ᕪ ᕫ ᕬ ᕭ ᕮ ᕯ ᕰ ᕱ ᕲ ᕳ ᕴ ᕵ ᕶ ᕷ ᕸ ᕹ ᕺ ᕻ ᕼ ᕽ ᕾ ᕿ
ᖀ ᖁ ᖂ ᖃ ᖄ ᖅ ᖆ ᖇ ᖈ ᖉ ᖊ ᖋ ᖌ ᖍ ᖎ ᖏ ᖐ ᖑ ᖒ ᖓ ᖔ ᖕ ᖖ ᖗ ᖘ ᖙ ᖚ ᖛ ᖜ ᖝ ᖞ ᖟ
ᖠ ᖡ ᖢ ᖣ ᖤ ᖥ ᖦ ᖧ ᖨ ᖩ ᖪ ᖫ ᖬ ᖭ ᖮ ᖯ ᖰ ᖱ ᖲ ᖳ ᖴ ᖵ ᖶ ᖷ ᖸ ᖹ ᖺ ᖻ ᖼ ᖽ ᖾ ᖿ
ᗀ ᗁ ᗂ ᗃ ᗄ ᗅ ᗆ ᗇ ᗈ ᗉ ᗊ ᗋ ᗌ ᗍ ᗎ ᗏ ᗐ ᗑ ᗒ ᗓ ᗔ ᗕ ᗖ ᗗ ᗘ ᗙ ᗚ ᗛ ᗜ ᗝ ᗞ ᗟ
ᗠ ᗡ ᗢ ᗣ ᗤ ᗥ ᗦ ᗧ ᗨ ᗩ ᗪ ᗫ ᗬ ᗭ ᗮ ᗯ ᗰ ᗱ ᗲ ᗳ ᗴ ᗵ ᗶ ᗷ ᗸ ᗹ ᗺ ᗻ ᗼ ᗽ ᗾ ᗿ
ᘀ ᘁ ᘂ ᘃ ᘄ ᘅ ᘆ ᘇ ᘈ ᘉ ᘊ ᘋ ᘌ ᘍ ᘎ ᘏ ᘐ ᘑ ᘒ ᘓ ᘔ ᘕ ᘖ ᘗ ᘘ ᘙ ᘚ ᘛ ᘜ ᘝ ᘞ ᘟ
ᘠ ᘡ ᘢ ᘣ ᘤ ᘥ ᘦ ᘧ ᘨ ᘩ ᘪ ᘫ ᘬ ᘭ ᘮ ᘯ ᘰ ᘱ ᘲ ᘳ ᘴ ᘵ ᘶ ᘷ ᘸ ᘹ ᘺ ᘻ ᘼ ᘽ ᘾ ᘿ
ᙀ ᙁ ᙂ ᙃ ᙄ ᙅ ᙆ ᙇ ᙈ ᙉ ᙊ ᙋ ᙌ ᙍ ᙎ ᙏ ᙐ ᙑ ᙒ ᙓ ᙔ ᙕ ᙖ ᙗ ᙘ ᙙ ᙚ ᙛ ᙜ ᙝ ᙞ ᙟ
ᙠ ᙡ ᙢ ᙣ ᙤ ᙥ ᙦ ᙧ ᙨ ᙩ ᙪ ᙫ ᙬ ᙭ ᙮ ᙯ ᙰ ᙱ ᙲ ᙳ ᙴ ᙵ ᙶ ᙷ ᙸ ᙹ ᙺ ᙻ ᙼ ᙽ ᙾ ᙿ
ᢰ ᢱ ᢲ ᢳ ᢴ ᢵ ᢶ ᢷ ᢸ ᢹ ᢺ ᢻ ᢼ ᢽ ᢾ ᢿ ᣀ ᣁ ᣂ ᣃ ᣄ ᣅ ᣆ ᣇ ᣈ ᣉ ᣊ ᣋ ᣌ ᣍ ᣎ ᣏ
ᣐ ᣑ ᣒ ᣓ ᣔ ᣕ ᣖ ᣗ ᣘ ᣙ ᣚ ᣛ ᣜ ᣝ ᣞ ᣟ ᣠ ᣡ ᣢ ᣣ ᣤ ᣥ ᣦ ᣧ ᣨ ᣩ ᣪ ᣫ ᣬ ᣭ ᣮ ᣯ
ᣰ ᣱ ᣲ ᣳ ᣴ ᣵ
""".strip()


def get_font(size=24):
    """Try to find a suitable font for Canadian Aboriginal Syllabics."""
    # Try common fonts with good Unicode support
    font_options = [
        "/usr/share/fonts/truetype/freefont/FreeSans.ttf",
        "DroidSansMono.ttf",
        "DejaVuSansMono.ttf",
        "NotoSans-Regular.ttf",
        "FreeMono.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf",
        "/usr/share/fonts/truetype/droid/DroidSansMono.ttf",
        "/usr/share/fonts/truetype/liberation/LiberationMono-Regular.ttf",
        "/usr/share/fonts/truetype/noto/NotoSans-Regular.ttf",
    ]

    for font_path in font_options:
        try:
            font = ImageFont.truetype(font_path, size)
            return font
        except (OSError, IOError):
            continue

    # Fallback to default font
    print("Warning: Could not find preferred font, using default")
    return ImageFont.load_default()

def screen_fonts(size = 24):
    
    test_char = "ᙗ"
    
    fontdir = Path("/usr/share/fonts/truetype")
    
    for font_fam in fontdir.iterdir():
        
        if not font_fam.is_dir():
            continue
        
        for font_file in font_fam.iterdir():
            
            if font_file.suffix != ".ttf":
                continue
            
            try:
                font = ImageFont.truetype(font_file, size = 32)
            except:
                print(f"failed to load font from {font_file.name}")
                continue
            
            rendered = render_character(test_char, size = 64, font = font)
            
            img_fname = f"./fonts/{font_file.stem}_test.png"            
            rendered.save(img_fname)


def render_character(char, size=48, font_size=36, font = None):
    """Render a single character as a PNG image."""
    # Create white background image
    img = Image.new('RGB', (size, size), color='white')
    draw = ImageDraw.Draw(img)

    if not font:
        # Load font
        font = get_font(font_size)

    # Get character bounding box for centering
    bbox = draw.textbbox((0, 0), char, font=font)
    text_width = bbox[2] - bbox[0]
    text_height = bbox[3] - bbox[1]

    # Center the character
    x = (size - text_width) // 2 - bbox[0]
    y = (size - text_height) // 2 - bbox[1]

    # Draw character in black
    draw.text((x, y), char, fill='black', font=font)

    return img


def render_all_characters(output_dir, size=48):
    """Render all characters and save them as individual PNGs."""
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True, parents=True)

    # Parse characters
    chars = CHARACTERS.split()

    print(f"Rendering {len(chars)} characters to {output_dir}/")

    for char in chars:
        # Get Unicode codepoint
        codepoint = ord(char)
        hex_code = f"U+{codepoint:04X}"

        # Render character
        img = render_character(char, size=size)

        # Save with descriptive filename
        # filename = f"{hex_code}_{char}.png"
        filename = f"{hex_code}.png"
        filepath = output_path / filename
        img.save(filepath)

        # print(f"  {hex_code} {char} -> {filename}")

    print(f"\nDone! Rendered {len(chars)} characters.")
    return chars


def create_grid_view(chars, output_dir, cols=16, cell_size=48):
    """Create a grid view showing all characters for easy reference."""
    rows = (len(chars) + cols - 1) // cols

    # Create large image for grid
    grid_width = cols * cell_size
    grid_height = rows * cell_size
    grid_img = Image.new('RGB', (grid_width, grid_height), color='white')

    font = get_font(36)

    for idx, char in enumerate(chars):
        row = idx // cols
        col = idx % cols

        # Render character
        char_img = render_character(char, size=cell_size, font_size=36)
        
        # Paste into grid
        x = col * cell_size
        y = row * cell_size
        grid_img.paste(char_img, (x, y))

    # Save grid
    grid_path = Path(output_dir) / "grid_view.png"
    grid_img.save(grid_path)
    print(f"\nGrid view saved to {grid_path}")

    return grid_img



def main():
    """Main function to render all characters."""
    
    # Output directory
    output_dir = "syllabics_rendered"

    # Render individual characters
    chars = render_all_characters(output_dir, size=48)

    # Create grid view
    create_grid_view(chars, output_dir, cols=32, cell_size=48)

    print(f"\n✓ All done! Check {output_dir}/ for rendered characters.")


if __name__ == "__main__":
    main()
