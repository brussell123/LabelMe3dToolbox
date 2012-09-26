function mesh = LM3DgetTexturedMesh(annotation,img)

mesh = getSceneMesh(annotation);
mesh = subDivideMesh(mesh,2);
mesh = getTextures(mesh,annotation,img);
