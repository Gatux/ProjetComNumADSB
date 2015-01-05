function [  ] = plane_pos_display( registre )


    PLANE_LON = registre.longitude % Longitude de l'avion
    PLANE_LAT = registre.latitude % Latitude de l'avion

    Id_airplane = 'plane';
    
    if isempty(registre.nom)
        Id_airplane = registre.adresse; % Nom de l'avion
    else
        Id_airplane = registre.nom;
    end
        
    if ~isempty(PLANE_LON) && ~isempty(PLANE_LAT)
        hold on
        plot(PLANE_LON,PLANE_LAT,'+b','MarkerSize',8);
        text(PLANE_LON+0.05,PLANE_LAT, num2str(Id_airplane),'color','b');
    end
    
    if ~isempty(registre.trajectoire)
        line([registre.trajectoire(1) PLANE_LON], [registre.trajectoire(2) PLANE_LAT]);
    end
    
end