const filename = "SupplementalTable_2_subset.csv";
const blob = new Blob([source], { type: "text/csv;charset=utf-8;" });

//addresses IE
if (navigator.msSaveBlob) {
  navigator.msSaveBlob(blob, filename);
} else {
  const link = document.createElement("a");
  link.href = URL.createObjectURL(blob);
  link.download = filename;
  link.target = "_blank";
  link.style.visibility = "hidden";
  link.dispatchEvent(new MouseEvent("click"));
}
