function make_video(working_dir, video_name, opt)
    % Generates a video starting from the files contained in the path
    % 'working_dir'.

    %% default options
    opt.dummy = 0;
    if ~isfield(opt,'frame_rate')
        opt.frame_rate = 25; % [fps]
    end
    
    %% make video
    outputVideo = VideoWriter(video_name);
    outputVideo.FrameRate = opt.frame_rate;
    %outputVideo.Quality = 50;
    open(outputVideo)

    imageNames = dir(fullfile(working_dir,'*.png'));
    imageNames = {imageNames.name}';
    for ii = 1:length(imageNames)
       fprintf('Image %d of %d...\n',ii,length(imageNames));
       img = imread(fullfile(working_dir,imageNames{ii}));
       writeVideo(outputVideo,img)
    end

    close(outputVideo)
    fprintf('video created!\n')
end