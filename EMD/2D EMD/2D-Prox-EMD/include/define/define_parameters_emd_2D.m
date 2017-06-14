function lambda = define_parameters_emd_2D(example,type)

if strcmp(example,'example1')
    if strcmp(type,'G2D')
        lambda.lambda_geometry{1} = 0.9;
        lambda.lambda_texture{1} = 1000;
    elseif strcmp(type,'G2Ddir')
        lambda.lambda_geometry{1} = 1.1;
        lambda.lambda_texture{1} = 1500;
    elseif strcmp(type,'P2DLC')
        lambda.lambda_geometry{1} = 2;
        lambda.lambda_texture_l{1} = 1;
        lambda.lambda_texture_c{1} = 1;
    elseif strcmp(type,'P2DLCD')
        lambda.lambda_geometry{1} = 2;
        lambda.lambda_texture_l{1} = 1;
        lambda.lambda_texture_c{1} = 1;
        lambda.lambda_texture_d1{1} = 1;
        lambda.lambda_texture_d2{1} = 1;
    end
   lambda.K=1;     
end
if strcmp(example,'example2')
    if strcmp(type,'G2D')
        lambda.lambda_geometry{1} = 0.1;
        lambda.lambda_texture{1} = 5;
        lambda.lambda_geometry{2} = 0.1;
        lambda.lambda_texture{2} = 0.001;
    elseif strcmp(type,'G2Ddir')
        lambda.lambda_geometry{1} = 0.1;
        lambda.lambda_texture{1} = 10;
        lambda.lambda_geometry{2} = 1;
        lambda.lambda_texture{2} = 0.1;
    elseif strcmp(type,'P2DLC')
        lambda.lambda_geometry{1} = 1;
        lambda.lambda_texture_l{1} = 1;
        lambda.lambda_texture_c{1} = 1;
        lambda.lambda_geometry{2} = 6;
        lambda.lambda_texture_l{2} = 1;
        lambda.lambda_texture_c{2} = 1;
    elseif strcmp(type,'P2DLCD')
        lambda.lambda_geometry{1} = 1;
        lambda.lambda_texture_l{1} = 1;
        lambda.lambda_texture_c{1} = 1;
        lambda.lambda_texture_d1{1} = 1;
        lambda.lambda_texture_d2{1} = 1;
        lambda.lambda_geometry{2} = 21;

        lambda.lambda_texture_l{2} = 1;
        lambda.lambda_texture_c{2} = 1;
        lambda.lambda_texture_d1{2} = 1;
        lambda.lambda_texture_d2{2} = 1;
    end
   lambda.K=2;     
end
if strcmp(example,'example3')
    if strcmp(type,'G2D')
        lambda.lambda_geometry{1} = 0.5;
        lambda.lambda_texture{1} = 100;
        lambda.lambda_geometry{2} = 0.05;
        lambda.lambda_texture{2} = 1;
    elseif strcmp(type,'G2Ddir')
        lambda.lambda_geometry{1} = 0.1;
        lambda.lambda_texture{1} = 200;
        lambda.lambda_geometry{2} = 0.05;
        lambda.lambda_texture{2} = 1;
    elseif strcmp(type,'P2DLC')
        lambda.lambda_geometry{1} = 0.3;
        lambda.lambda_texture_l{1} = 0.5;
        lambda.lambda_texture_c{1} = 0.5;
        lambda.lambda_geometry{2} = 1;
        lambda.lambda_texture_l{2} = 0.1;
        lambda.lambda_texture_c{2} = 0.1;
    elseif strcmp(type,'P2DLCD')
        lambda.lambda_geometry{1} = 0.3;
        lambda.lambda_texture_l{1} = 0.3;
        lambda.lambda_texture_c{1} = 0.3;
        lambda.lambda_texture_d1{1} = 0.3;
        lambda.lambda_texture_d2{1} = 0.3;
        lambda.lambda_geometry{2} = 1;
        lambda.lambda_texture_l{2} = 0.1;
        lambda.lambda_texture_c{2} = 0.1;
        lambda.lambda_texture_d1{2} = 0.1;
        lambda.lambda_texture_d2{2} = 0.1;
    end
   lambda.K=2;    
end

end
