import kipoi
import kipoi.model
import kipoi.pipeline

# model = kipoi.get_model("DeepBind/Homo_sapiens/TF/D00328.018_ChIP-seq_CTCF")
model = kipoi.get_model("Basset")
# model = kipoi.get_model("Basset", source='kipoi', with_dataloader=True)
# pred = model.pipeline.predict_example(batch_size=4)