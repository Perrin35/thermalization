OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6948833) q[0];
sx q[0];
rz(2.0599685) q[0];
sx q[0];
rz(10.164227) q[0];
rz(2.3820355) q[1];
sx q[1];
rz(-1.3146725) q[1];
sx q[1];
rz(1.4256328) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.254296) q[0];
sx q[0];
rz(-2.097541) q[0];
sx q[0];
rz(0.33381427) q[0];
rz(2.6447634) q[2];
sx q[2];
rz(-1.3303927) q[2];
sx q[2];
rz(-0.59076004) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.4233746) q[1];
sx q[1];
rz(-1.6483232) q[1];
sx q[1];
rz(1.7192057) q[1];
rz(-pi) q[2];
rz(1.5208779) q[3];
sx q[3];
rz(-0.71943362) q[3];
sx q[3];
rz(-2.0572544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.90929675) q[2];
sx q[2];
rz(-2.2984419) q[2];
sx q[2];
rz(2.9006531) q[2];
rz(-3.1124034) q[3];
sx q[3];
rz(-1.3386644) q[3];
sx q[3];
rz(1.9624814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9772684) q[0];
sx q[0];
rz(-1.203953) q[0];
sx q[0];
rz(-0.48467317) q[0];
rz(-2.4618705) q[1];
sx q[1];
rz(-1.2815963) q[1];
sx q[1];
rz(-1.1601123) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74014464) q[0];
sx q[0];
rz(-0.30456802) q[0];
sx q[0];
rz(-0.061621678) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4258575) q[2];
sx q[2];
rz(-2.1308427) q[2];
sx q[2];
rz(-0.83362388) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.94494941) q[1];
sx q[1];
rz(-1.8183876) q[1];
sx q[1];
rz(-0.61243122) q[1];
rz(-pi) q[2];
rz(1.9002565) q[3];
sx q[3];
rz(-1.0715535) q[3];
sx q[3];
rz(-1.3105621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.844187) q[2];
sx q[2];
rz(-2.83941) q[2];
sx q[2];
rz(-0.45787946) q[2];
rz(1.1229905) q[3];
sx q[3];
rz(-2.1419958) q[3];
sx q[3];
rz(-0.56639731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26894012) q[0];
sx q[0];
rz(-2.432423) q[0];
sx q[0];
rz(3.1100682) q[0];
rz(0.28745502) q[1];
sx q[1];
rz(-0.87702409) q[1];
sx q[1];
rz(-1.215975) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8376802) q[0];
sx q[0];
rz(-1.0420351) q[0];
sx q[0];
rz(-0.75334511) q[0];
rz(0.38162614) q[2];
sx q[2];
rz(-0.51593057) q[2];
sx q[2];
rz(-0.53781539) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.79920125) q[1];
sx q[1];
rz(-1.2974129) q[1];
sx q[1];
rz(2.8552516) q[1];
x q[2];
rz(0.76283703) q[3];
sx q[3];
rz(-2.3194072) q[3];
sx q[3];
rz(-1.959182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.03881255) q[2];
sx q[2];
rz(-1.5993885) q[2];
sx q[2];
rz(2.3056324) q[2];
rz(-1.3577667) q[3];
sx q[3];
rz(-2.416555) q[3];
sx q[3];
rz(-2.3939705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8589856) q[0];
sx q[0];
rz(-1.7361807) q[0];
sx q[0];
rz(-3.0614241) q[0];
rz(-1.0824341) q[1];
sx q[1];
rz(-0.23324649) q[1];
sx q[1];
rz(-1.3287883) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4565312) q[0];
sx q[0];
rz(-1.0624806) q[0];
sx q[0];
rz(-1.1969181) q[0];
rz(-pi) q[1];
rz(-2.8187739) q[2];
sx q[2];
rz(-1.4070961) q[2];
sx q[2];
rz(1.1522918) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0094686) q[1];
sx q[1];
rz(-1.7156202) q[1];
sx q[1];
rz(-2.0028466) q[1];
x q[2];
rz(1.7747709) q[3];
sx q[3];
rz(-1.9001203) q[3];
sx q[3];
rz(2.9993111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3749915) q[2];
sx q[2];
rz(-2.1583755) q[2];
sx q[2];
rz(-0.90744606) q[2];
rz(0.67029101) q[3];
sx q[3];
rz(-2.3136316) q[3];
sx q[3];
rz(-1.420174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9012673) q[0];
sx q[0];
rz(-2.2690052) q[0];
sx q[0];
rz(0.43148828) q[0];
rz(-1.1314499) q[1];
sx q[1];
rz(-1.9624458) q[1];
sx q[1];
rz(-0.44050899) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90135306) q[0];
sx q[0];
rz(-0.47751891) q[0];
sx q[0];
rz(1.8285455) q[0];
rz(-2.0149286) q[2];
sx q[2];
rz(-2.6878549) q[2];
sx q[2];
rz(-2.1778088) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2825644) q[1];
sx q[1];
rz(-0.33278123) q[1];
sx q[1];
rz(3.0206693) q[1];
rz(-pi) q[2];
x q[2];
rz(0.93705658) q[3];
sx q[3];
rz(-1.0770633) q[3];
sx q[3];
rz(0.99668324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.12209192) q[2];
sx q[2];
rz(-2.2152405) q[2];
sx q[2];
rz(-2.7276373) q[2];
rz(-1.3881989) q[3];
sx q[3];
rz(-1.5475169) q[3];
sx q[3];
rz(0.12510124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6269161) q[0];
sx q[0];
rz(-1.4056982) q[0];
sx q[0];
rz(-0.93389121) q[0];
rz(-2.4712708) q[1];
sx q[1];
rz(-1.8972081) q[1];
sx q[1];
rz(1.3353039) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40437296) q[0];
sx q[0];
rz(-1.4708733) q[0];
sx q[0];
rz(2.1735682) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9839758) q[2];
sx q[2];
rz(-2.981957) q[2];
sx q[2];
rz(0.0052513382) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.28999968) q[1];
sx q[1];
rz(-1.4581175) q[1];
sx q[1];
rz(3.1188909) q[1];
rz(-pi) q[2];
rz(-1.5794831) q[3];
sx q[3];
rz(-0.93378769) q[3];
sx q[3];
rz(0.52906936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2495217) q[2];
sx q[2];
rz(-2.8688909) q[2];
sx q[2];
rz(-1.8708694) q[2];
rz(1.3696085) q[3];
sx q[3];
rz(-1.9000051) q[3];
sx q[3];
rz(0.52186596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9524566) q[0];
sx q[0];
rz(-0.58681762) q[0];
sx q[0];
rz(0.15368803) q[0];
rz(0.20248374) q[1];
sx q[1];
rz(-1.9053562) q[1];
sx q[1];
rz(0.94857803) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62371333) q[0];
sx q[0];
rz(-0.9615295) q[0];
sx q[0];
rz(-1.3246956) q[0];
rz(0.3237299) q[2];
sx q[2];
rz(-2.1846131) q[2];
sx q[2];
rz(-0.15440369) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7708009) q[1];
sx q[1];
rz(-1.8026514) q[1];
sx q[1];
rz(1.8370017) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7465215) q[3];
sx q[3];
rz(-2.4843289) q[3];
sx q[3];
rz(-2.6906309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1367246) q[2];
sx q[2];
rz(-2.0437045) q[2];
sx q[2];
rz(3.1401805) q[2];
rz(-3.0794365) q[3];
sx q[3];
rz(-1.4096189) q[3];
sx q[3];
rz(0.96127659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.38178) q[0];
sx q[0];
rz(-0.75730046) q[0];
sx q[0];
rz(-1.9267474) q[0];
rz(-0.75792056) q[1];
sx q[1];
rz(-2.6140723) q[1];
sx q[1];
rz(-2.5835999) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.081462534) q[0];
sx q[0];
rz(-1.234982) q[0];
sx q[0];
rz(2.8804638) q[0];
rz(-pi) q[1];
rz(-2.5900389) q[2];
sx q[2];
rz(-1.4467014) q[2];
sx q[2];
rz(-2.6089904) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3000189) q[1];
sx q[1];
rz(-1.2928559) q[1];
sx q[1];
rz(-2.071408) q[1];
rz(-1.6648983) q[3];
sx q[3];
rz(-2.8633139) q[3];
sx q[3];
rz(-0.32839963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.5672292) q[2];
sx q[2];
rz(-1.9476674) q[2];
sx q[2];
rz(2.8723259) q[2];
rz(-1.1632129) q[3];
sx q[3];
rz(-2.5389157) q[3];
sx q[3];
rz(2.9009624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4686541) q[0];
sx q[0];
rz(-2.1042295) q[0];
sx q[0];
rz(1.0733676) q[0];
rz(2.0138373) q[1];
sx q[1];
rz(-2.914371) q[1];
sx q[1];
rz(-2.5849297) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8506691) q[0];
sx q[0];
rz(-1.2525096) q[0];
sx q[0];
rz(-0.14936635) q[0];
x q[1];
rz(-2.3586541) q[2];
sx q[2];
rz(-0.73854827) q[2];
sx q[2];
rz(-0.99817456) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.6491739) q[1];
sx q[1];
rz(-0.88090501) q[1];
sx q[1];
rz(1.0997195) q[1];
rz(-0.39770522) q[3];
sx q[3];
rz(-0.25358646) q[3];
sx q[3];
rz(1.2802779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2399981) q[2];
sx q[2];
rz(-2.0549213) q[2];
sx q[2];
rz(1.1617917) q[2];
rz(2.6084172) q[3];
sx q[3];
rz(-0.52459255) q[3];
sx q[3];
rz(-0.1575135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48851442) q[0];
sx q[0];
rz(-2.1670659) q[0];
sx q[0];
rz(0.59481204) q[0];
rz(-0.57890233) q[1];
sx q[1];
rz(-1.4756823) q[1];
sx q[1];
rz(-0.65931177) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66856495) q[0];
sx q[0];
rz(-0.94946948) q[0];
sx q[0];
rz(-1.4061808) q[0];
rz(2.012678) q[2];
sx q[2];
rz(-1.0143447) q[2];
sx q[2];
rz(2.93612) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.54592268) q[1];
sx q[1];
rz(-2.8575142) q[1];
sx q[1];
rz(2.0998663) q[1];
x q[2];
rz(1.8439383) q[3];
sx q[3];
rz(-1.3160664) q[3];
sx q[3];
rz(-1.8984853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.79002964) q[2];
sx q[2];
rz(-1.943925) q[2];
sx q[2];
rz(-0.053038049) q[2];
rz(1.3430345) q[3];
sx q[3];
rz(-1.2767295) q[3];
sx q[3];
rz(0.0742577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85528436) q[0];
sx q[0];
rz(-1.3636148) q[0];
sx q[0];
rz(-1.6028945) q[0];
rz(3.042649) q[1];
sx q[1];
rz(-1.7715441) q[1];
sx q[1];
rz(-3.0861707) q[1];
rz(2.7361353) q[2];
sx q[2];
rz(-1.5787072) q[2];
sx q[2];
rz(-2.2070259) q[2];
rz(-2.4002038) q[3];
sx q[3];
rz(-1.6869873) q[3];
sx q[3];
rz(1.6968873) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
