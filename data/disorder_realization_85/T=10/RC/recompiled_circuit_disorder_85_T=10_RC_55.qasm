OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.20733362) q[0];
sx q[0];
rz(-2.5512295) q[0];
sx q[0];
rz(-0.37101775) q[0];
rz(-0.38129216) q[1];
sx q[1];
rz(-0.59950221) q[1];
sx q[1];
rz(-1.7655656) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3289514) q[0];
sx q[0];
rz(-2.0094123) q[0];
sx q[0];
rz(0.89226512) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0618093) q[2];
sx q[2];
rz(-2.2248189) q[2];
sx q[2];
rz(2.430928) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.10378097) q[1];
sx q[1];
rz(-1.4191287) q[1];
sx q[1];
rz(1.9019466) q[1];
x q[2];
rz(2.3222378) q[3];
sx q[3];
rz(-1.5324394) q[3];
sx q[3];
rz(-1.0591782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.084289) q[2];
sx q[2];
rz(-0.40033445) q[2];
sx q[2];
rz(-0.98891813) q[2];
rz(0.75254285) q[3];
sx q[3];
rz(-1.1458594) q[3];
sx q[3];
rz(-0.83077103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20724021) q[0];
sx q[0];
rz(-3.0293284) q[0];
sx q[0];
rz(1.1799312) q[0];
rz(2.143899) q[1];
sx q[1];
rz(-1.2832063) q[1];
sx q[1];
rz(2.4172799) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1441919) q[0];
sx q[0];
rz(-1.8840944) q[0];
sx q[0];
rz(-2.2749167) q[0];
x q[1];
rz(1.486859) q[2];
sx q[2];
rz(-1.3720023) q[2];
sx q[2];
rz(2.4059911) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.95485479) q[1];
sx q[1];
rz(-1.2192982) q[1];
sx q[1];
rz(-0.21912205) q[1];
x q[2];
rz(2.4507387) q[3];
sx q[3];
rz(-2.8375531) q[3];
sx q[3];
rz(-0.090304852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.90536845) q[2];
sx q[2];
rz(-1.3304109) q[2];
sx q[2];
rz(0.36188564) q[2];
rz(0.13606717) q[3];
sx q[3];
rz(-2.5858904) q[3];
sx q[3];
rz(0.045624174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9996465) q[0];
sx q[0];
rz(-1.0936341) q[0];
sx q[0];
rz(1.746159) q[0];
rz(-2.6793001) q[1];
sx q[1];
rz(-0.42458436) q[1];
sx q[1];
rz(1.9225072) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.483706) q[0];
sx q[0];
rz(-1.833796) q[0];
sx q[0];
rz(-1.2067243) q[0];
rz(-2.8457882) q[2];
sx q[2];
rz(-0.89106262) q[2];
sx q[2];
rz(0.55065292) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3238941) q[1];
sx q[1];
rz(-2.5496799) q[1];
sx q[1];
rz(1.1281668) q[1];
rz(-pi) q[2];
x q[2];
rz(1.838802) q[3];
sx q[3];
rz(-1.5857113) q[3];
sx q[3];
rz(-2.5915495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7200155) q[2];
sx q[2];
rz(-0.74024671) q[2];
sx q[2];
rz(1.4455618) q[2];
rz(-0.56882632) q[3];
sx q[3];
rz(-2.2918662) q[3];
sx q[3];
rz(-0.11051699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22359426) q[0];
sx q[0];
rz(-0.44819865) q[0];
sx q[0];
rz(-0.51825994) q[0];
rz(-2.4261684) q[1];
sx q[1];
rz(-2.0253069) q[1];
sx q[1];
rz(-2.3148361) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0543538) q[0];
sx q[0];
rz(-1.0687437) q[0];
sx q[0];
rz(-1.0994083) q[0];
x q[1];
rz(-1.4364169) q[2];
sx q[2];
rz(-2.3339286) q[2];
sx q[2];
rz(-0.32854167) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3863694) q[1];
sx q[1];
rz(-1.1896903) q[1];
sx q[1];
rz(-1.2632881) q[1];
rz(-pi) q[2];
rz(0.34927807) q[3];
sx q[3];
rz(-0.75540245) q[3];
sx q[3];
rz(-0.7790156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9892019) q[2];
sx q[2];
rz(-2.9426136) q[2];
sx q[2];
rz(-1.7626804) q[2];
rz(0.072323024) q[3];
sx q[3];
rz(-0.81243378) q[3];
sx q[3];
rz(-1.6453843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3146661) q[0];
sx q[0];
rz(-2.5214654) q[0];
sx q[0];
rz(2.0157053) q[0];
rz(0.90244883) q[1];
sx q[1];
rz(-2.1676962) q[1];
sx q[1];
rz(-2.856423) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0774768) q[0];
sx q[0];
rz(-2.6649094) q[0];
sx q[0];
rz(1.478273) q[0];
x q[1];
rz(-0.4816149) q[2];
sx q[2];
rz(-2.4184605) q[2];
sx q[2];
rz(2.8358104) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.8385182) q[1];
sx q[1];
rz(-2.8285366) q[1];
sx q[1];
rz(1.5721553) q[1];
rz(1.6676335) q[3];
sx q[3];
rz(-1.9707465) q[3];
sx q[3];
rz(-1.4777583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.29331648) q[2];
sx q[2];
rz(-0.50555503) q[2];
sx q[2];
rz(2.6021393) q[2];
rz(-2.8347677) q[3];
sx q[3];
rz(-2.2570733) q[3];
sx q[3];
rz(2.6873798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8834615) q[0];
sx q[0];
rz(-2.3911609) q[0];
sx q[0];
rz(-0.15701292) q[0];
rz(2.4482588) q[1];
sx q[1];
rz(-0.88070977) q[1];
sx q[1];
rz(1.7745811) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0317504) q[0];
sx q[0];
rz(-1.5378008) q[0];
sx q[0];
rz(1.6389636) q[0];
rz(-pi) q[1];
x q[1];
rz(0.579367) q[2];
sx q[2];
rz(-0.89951347) q[2];
sx q[2];
rz(-2.6210149) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4821266) q[1];
sx q[1];
rz(-0.30059338) q[1];
sx q[1];
rz(1.8468922) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.35297024) q[3];
sx q[3];
rz(-1.4735231) q[3];
sx q[3];
rz(-0.73247611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.29629016) q[2];
sx q[2];
rz(-2.2915816) q[2];
sx q[2];
rz(-0.40346754) q[2];
rz(0.48163313) q[3];
sx q[3];
rz(-1.0721595) q[3];
sx q[3];
rz(-0.51923716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24213174) q[0];
sx q[0];
rz(-2.2583028) q[0];
sx q[0];
rz(-0.8738628) q[0];
rz(2.6938687) q[1];
sx q[1];
rz(-0.73900765) q[1];
sx q[1];
rz(1.1706932) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51629492) q[0];
sx q[0];
rz(-1.5769616) q[0];
sx q[0];
rz(-0.023029285) q[0];
rz(-pi) q[1];
rz(0.41202338) q[2];
sx q[2];
rz(-0.12672666) q[2];
sx q[2];
rz(0.7574946) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.891174) q[1];
sx q[1];
rz(-2.182057) q[1];
sx q[1];
rz(2.3120018) q[1];
rz(-pi) q[2];
rz(0.4865173) q[3];
sx q[3];
rz(-0.58972893) q[3];
sx q[3];
rz(-0.041134838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0968904) q[2];
sx q[2];
rz(-2.5723852) q[2];
sx q[2];
rz(0.61075413) q[2];
rz(-2.6664873) q[3];
sx q[3];
rz(-2.0510309) q[3];
sx q[3];
rz(0.92774123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24191813) q[0];
sx q[0];
rz(-3.0245259) q[0];
sx q[0];
rz(-0.29712594) q[0];
rz(1.3946474) q[1];
sx q[1];
rz(-1.9906094) q[1];
sx q[1];
rz(-0.64613211) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.811759) q[0];
sx q[0];
rz(-0.23010294) q[0];
sx q[0];
rz(2.0287201) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8201639) q[2];
sx q[2];
rz(-0.49389631) q[2];
sx q[2];
rz(2.0733881) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7039973) q[1];
sx q[1];
rz(-0.62347177) q[1];
sx q[1];
rz(2.5295579) q[1];
rz(-pi) q[2];
x q[2];
rz(2.331278) q[3];
sx q[3];
rz(-1.1018361) q[3];
sx q[3];
rz(-1.9988434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0769161) q[2];
sx q[2];
rz(-0.95305324) q[2];
sx q[2];
rz(2.6867552) q[2];
rz(0.70139766) q[3];
sx q[3];
rz(-2.111179) q[3];
sx q[3];
rz(2.0075683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4421473) q[0];
sx q[0];
rz(-3*pi/16) q[0];
sx q[0];
rz(-0.79750693) q[0];
rz(0.51756716) q[1];
sx q[1];
rz(-2.3289754) q[1];
sx q[1];
rz(0.10841766) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3603044) q[0];
sx q[0];
rz(-0.66626781) q[0];
sx q[0];
rz(-1.4772619) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3181395) q[2];
sx q[2];
rz(-2.9571819) q[2];
sx q[2];
rz(2.2968963) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.11549599) q[1];
sx q[1];
rz(-2.8671088) q[1];
sx q[1];
rz(-0.91067578) q[1];
rz(-pi) q[2];
rz(2.9547144) q[3];
sx q[3];
rz(-2.2664321) q[3];
sx q[3];
rz(2.7587492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1291528) q[2];
sx q[2];
rz(-1.3468578) q[2];
sx q[2];
rz(-0.33995315) q[2];
rz(2.7231976) q[3];
sx q[3];
rz(-0.59643006) q[3];
sx q[3];
rz(-2.4160014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6091992) q[0];
sx q[0];
rz(-0.39390716) q[0];
sx q[0];
rz(0.6788196) q[0];
rz(2.7774096) q[1];
sx q[1];
rz(-1.4441676) q[1];
sx q[1];
rz(3.0864339) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63294166) q[0];
sx q[0];
rz(-1.5061437) q[0];
sx q[0];
rz(2.7642194) q[0];
rz(0.55969413) q[2];
sx q[2];
rz(-1.7494546) q[2];
sx q[2];
rz(-0.44911227) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0835417) q[1];
sx q[1];
rz(-1.0293048) q[1];
sx q[1];
rz(-0.70152775) q[1];
x q[2];
rz(2.0268029) q[3];
sx q[3];
rz(-1.5378693) q[3];
sx q[3];
rz(-1.959521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1577592) q[2];
sx q[2];
rz(-1.021421) q[2];
sx q[2];
rz(2.6514163) q[2];
rz(-0.13752078) q[3];
sx q[3];
rz(-1.0995882) q[3];
sx q[3];
rz(-2.2035051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72538439) q[0];
sx q[0];
rz(-1.2513456) q[0];
sx q[0];
rz(-0.67847897) q[0];
rz(2.9329119) q[1];
sx q[1];
rz(-1.122767) q[1];
sx q[1];
rz(-1.541419) q[1];
rz(0.31221496) q[2];
sx q[2];
rz(-1.4785462) q[2];
sx q[2];
rz(-1.2350456) q[2];
rz(-0.18110885) q[3];
sx q[3];
rz(-2.8977179) q[3];
sx q[3];
rz(0.44117622) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];