OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5118339) q[0];
sx q[0];
rz(-0.30183733) q[0];
sx q[0];
rz(-1.8344185) q[0];
rz(-3.6699927) q[1];
sx q[1];
rz(3.9730605) q[1];
sx q[1];
rz(12.107036) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49469264) q[0];
sx q[0];
rz(-1.6199058) q[0];
sx q[0];
rz(-1.3145334) q[0];
rz(-1.2940965) q[2];
sx q[2];
rz(-1.3802229) q[2];
sx q[2];
rz(2.0399567) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2832773) q[1];
sx q[1];
rz(-2.1766571) q[1];
sx q[1];
rz(-1.359904) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1189947) q[3];
sx q[3];
rz(-0.90211464) q[3];
sx q[3];
rz(2.150226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.022973148) q[2];
sx q[2];
rz(-1.8034673) q[2];
sx q[2];
rz(2.5845134) q[2];
rz(-2.6856375) q[3];
sx q[3];
rz(-2.3796701) q[3];
sx q[3];
rz(-2.5383811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93422455) q[0];
sx q[0];
rz(-0.71894431) q[0];
sx q[0];
rz(1.7508605) q[0];
rz(-1.8234183) q[1];
sx q[1];
rz(-1.7134106) q[1];
sx q[1];
rz(-2.9255829) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.073450449) q[0];
sx q[0];
rz(-1.8863858) q[0];
sx q[0];
rz(-0.95629779) q[0];
rz(-pi) q[1];
rz(-0.41925307) q[2];
sx q[2];
rz(-1.6800095) q[2];
sx q[2];
rz(-0.85604233) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.60571721) q[1];
sx q[1];
rz(-0.89369666) q[1];
sx q[1];
rz(0.89800055) q[1];
rz(-pi) q[2];
rz(2.9227867) q[3];
sx q[3];
rz(-2.1010733) q[3];
sx q[3];
rz(1.0007364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7184489) q[2];
sx q[2];
rz(-2.7687912) q[2];
sx q[2];
rz(0.049840363) q[2];
rz(0.54346624) q[3];
sx q[3];
rz(-1.1246357) q[3];
sx q[3];
rz(-1.0928833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8922888) q[0];
sx q[0];
rz(-2.0396621) q[0];
sx q[0];
rz(-2.1732543) q[0];
rz(-0.020542055) q[1];
sx q[1];
rz(-1.7612709) q[1];
sx q[1];
rz(1.4302018) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.641159) q[0];
sx q[0];
rz(-0.64167385) q[0];
sx q[0];
rz(3.1042751) q[0];
rz(-pi) q[1];
rz(1.9400575) q[2];
sx q[2];
rz(-1.0383237) q[2];
sx q[2];
rz(-2.7209501) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9544977) q[1];
sx q[1];
rz(-1.375838) q[1];
sx q[1];
rz(2.7379237) q[1];
x q[2];
rz(-2.7647721) q[3];
sx q[3];
rz(-1.6018036) q[3];
sx q[3];
rz(3.0135807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.7057544) q[2];
sx q[2];
rz(-0.64283723) q[2];
sx q[2];
rz(-1.359681) q[2];
rz(0.5082353) q[3];
sx q[3];
rz(-1.619092) q[3];
sx q[3];
rz(-1.1448418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0559167) q[0];
sx q[0];
rz(-2.8450232) q[0];
sx q[0];
rz(-1.0067518) q[0];
rz(2.309917) q[1];
sx q[1];
rz(-2.3995903) q[1];
sx q[1];
rz(-2.0337909) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9079893) q[0];
sx q[0];
rz(-0.74626669) q[0];
sx q[0];
rz(2.1265592) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8227656) q[2];
sx q[2];
rz(-1.2059847) q[2];
sx q[2];
rz(-0.081720086) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0084787) q[1];
sx q[1];
rz(-1.263887) q[1];
sx q[1];
rz(3.1208193) q[1];
x q[2];
rz(1.8516225) q[3];
sx q[3];
rz(-1.1410332) q[3];
sx q[3];
rz(-2.6971779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6882249) q[2];
sx q[2];
rz(-1.0658762) q[2];
sx q[2];
rz(0.31309703) q[2];
rz(-2.744216) q[3];
sx q[3];
rz(-1.671096) q[3];
sx q[3];
rz(-2.9562922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50936407) q[0];
sx q[0];
rz(-2.2780184) q[0];
sx q[0];
rz(-2.2160227) q[0];
rz(0.46982345) q[1];
sx q[1];
rz(-1.3146105) q[1];
sx q[1];
rz(2.125461) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6189656) q[0];
sx q[0];
rz(-2.5182808) q[0];
sx q[0];
rz(-0.81617426) q[0];
rz(1.2193331) q[2];
sx q[2];
rz(-1.8242852) q[2];
sx q[2];
rz(3.0712822) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.092246902) q[1];
sx q[1];
rz(-1.7723485) q[1];
sx q[1];
rz(1.3231897) q[1];
x q[2];
rz(-1.8941325) q[3];
sx q[3];
rz(-1.6010906) q[3];
sx q[3];
rz(2.7504675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4887345) q[2];
sx q[2];
rz(-2.8159339) q[2];
sx q[2];
rz(2.6540836) q[2];
rz(2.2440535) q[3];
sx q[3];
rz(-1.8833912) q[3];
sx q[3];
rz(-2.0770309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2010736) q[0];
sx q[0];
rz(-0.93795332) q[0];
sx q[0];
rz(2.4679389) q[0];
rz(1.2334088) q[1];
sx q[1];
rz(-1.0181095) q[1];
sx q[1];
rz(-0.76603755) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85114408) q[0];
sx q[0];
rz(-0.58935114) q[0];
sx q[0];
rz(-1.0714156) q[0];
x q[1];
rz(2.7577007) q[2];
sx q[2];
rz(-0.18927477) q[2];
sx q[2];
rz(1.4658907) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8173816) q[1];
sx q[1];
rz(-1.7158475) q[1];
sx q[1];
rz(-2.0562028) q[1];
rz(-pi) q[2];
rz(2.7286639) q[3];
sx q[3];
rz(-0.33187643) q[3];
sx q[3];
rz(1.6433034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.172714) q[2];
sx q[2];
rz(-1.4977027) q[2];
sx q[2];
rz(-0.76095757) q[2];
rz(2.739665) q[3];
sx q[3];
rz(-0.26508489) q[3];
sx q[3];
rz(-0.5886122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6589979) q[0];
sx q[0];
rz(-0.75356475) q[0];
sx q[0];
rz(-1.448451) q[0];
rz(-0.50085577) q[1];
sx q[1];
rz(-1.294699) q[1];
sx q[1];
rz(-0.81333152) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.017352176) q[0];
sx q[0];
rz(-2.725714) q[0];
sx q[0];
rz(0.37353596) q[0];
rz(-pi) q[1];
rz(-2.1189999) q[2];
sx q[2];
rz(-1.6728901) q[2];
sx q[2];
rz(0.36985794) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.2351027) q[1];
sx q[1];
rz(-2.0079327) q[1];
sx q[1];
rz(-1.4187995) q[1];
x q[2];
rz(1.4485095) q[3];
sx q[3];
rz(-1.5124195) q[3];
sx q[3];
rz(-0.98298847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8616051) q[2];
sx q[2];
rz(-0.6051175) q[2];
sx q[2];
rz(-2.81847) q[2];
rz(-1.7396287) q[3];
sx q[3];
rz(-1.2819382) q[3];
sx q[3];
rz(-1.3824979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.062227) q[0];
sx q[0];
rz(-1.0492188) q[0];
sx q[0];
rz(2.3754689) q[0];
rz(-0.8194204) q[1];
sx q[1];
rz(-2.1111646) q[1];
sx q[1];
rz(1.6105509) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19056828) q[0];
sx q[0];
rz(-1.6389209) q[0];
sx q[0];
rz(-0.22916746) q[0];
rz(-0.28240164) q[2];
sx q[2];
rz(-0.67495433) q[2];
sx q[2];
rz(-2.9815428) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.69680113) q[1];
sx q[1];
rz(-2.2292622) q[1];
sx q[1];
rz(-0.80689238) q[1];
rz(-1.2458477) q[3];
sx q[3];
rz(-1.2940931) q[3];
sx q[3];
rz(0.70223728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3766342) q[2];
sx q[2];
rz(-2.9464293) q[2];
sx q[2];
rz(2.7435319) q[2];
rz(-2.5352488) q[3];
sx q[3];
rz(-1.5747993) q[3];
sx q[3];
rz(2.2280367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23676087) q[0];
sx q[0];
rz(-2.8396711) q[0];
sx q[0];
rz(2.6171369) q[0];
rz(0.66871387) q[1];
sx q[1];
rz(-1.6122183) q[1];
sx q[1];
rz(-1.3800157) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3052914) q[0];
sx q[0];
rz(-2.143465) q[0];
sx q[0];
rz(3.089726) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0933268) q[2];
sx q[2];
rz(-1.8430955) q[2];
sx q[2];
rz(-1.7881025) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3375597) q[1];
sx q[1];
rz(-1.9849249) q[1];
sx q[1];
rz(-1.1131338) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5502167) q[3];
sx q[3];
rz(-2.5246361) q[3];
sx q[3];
rz(-2.5643831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1606719) q[2];
sx q[2];
rz(-2.394684) q[2];
sx q[2];
rz(-0.46372947) q[2];
rz(1.7175698) q[3];
sx q[3];
rz(-2.0537036) q[3];
sx q[3];
rz(3.1080642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56228191) q[0];
sx q[0];
rz(-2.4233241) q[0];
sx q[0];
rz(-0.41807362) q[0];
rz(-2.383291) q[1];
sx q[1];
rz(-0.98379358) q[1];
sx q[1];
rz(1.8565348) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83259634) q[0];
sx q[0];
rz(-0.65387883) q[0];
sx q[0];
rz(0.12419463) q[0];
rz(-1.1056771) q[2];
sx q[2];
rz(-0.67567247) q[2];
sx q[2];
rz(-3.1356406) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9662719) q[1];
sx q[1];
rz(-1.0283677) q[1];
sx q[1];
rz(1.7955417) q[1];
x q[2];
rz(2.140652) q[3];
sx q[3];
rz(-2.7469198) q[3];
sx q[3];
rz(-0.81854023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.16984223) q[2];
sx q[2];
rz(-1.0047793) q[2];
sx q[2];
rz(2.8954835) q[2];
rz(1.8698112) q[3];
sx q[3];
rz(-1.5747986) q[3];
sx q[3];
rz(-1.3625712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9825738) q[0];
sx q[0];
rz(-1.6652501) q[0];
sx q[0];
rz(-0.2022947) q[0];
rz(-0.62494878) q[1];
sx q[1];
rz(-2.2849871) q[1];
sx q[1];
rz(-2.7402592) q[1];
rz(1.3622147) q[2];
sx q[2];
rz(-1.7558395) q[2];
sx q[2];
rz(0.21370733) q[2];
rz(2.940964) q[3];
sx q[3];
rz(-1.8798141) q[3];
sx q[3];
rz(0.78165913) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
