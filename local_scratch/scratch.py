# Function to read the guide library file in .csv and .txt format
-


def read_guides_file(guides_file):
    -    file_extension = os.path.splitext(guides_file)[1].lower()


-
-
if file_extension == '.txt':
    -  # Read as a tab-delimited file
-        guides_df = pd.read_csv(guides_file, sep = '\t', engine = 'c')
- elif file_extension == '.csv':
-  # Read as a comma-separated file
-        guides_df = pd.read_csv(guides_file, engine = 'c')
- else:
-
raise ValueError("Unsupported file type: {}".format(file_extension))
+  # TODO: Ensure function runs
+  # def read_guides_file(guides_file):
+  # file_extension = os.path.splitext(guides_file)[1].lower()
+  #
+  # if file_extension == '.txt':
+  # # Read as a tab-delimited file
+  # guides_df = pd.read_csv(guides_file, sep = '\t', engine = 'c')
+  # elif file_extension == '.csv':
+  # # Read as a comma-separated file
+  # guides_df = pd.read_csv(guides_file, engine = 'c')
+  # else:
+  # raise ValueError("Unsupported file type: {}".format(file_extension))
+  #
+  # return guides_df
+  #
+  # guides_df = read_guides_file(guides_file)

-
return guides_df
-
-guides_df = read_guides_file(guides_file)
+  # Modify based on guide file format - Guide file received was in txt format
+guides_df = pd.read_csv(guides_file, sep = '\t')

"""# # Import R1s and R2s.
# ## pysam way Pysam introduces many other file format compatibilities such as
