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
rz(0.63737386) q[0];
sx q[0];
rz(-2.0847991) q[0];
sx q[0];
rz(-0.771653) q[0];
rz(-0.15089384) q[1];
sx q[1];
rz(-0.3124736) q[1];
sx q[1];
rz(-1.6627275) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2914198) q[0];
sx q[0];
rz(-1.2787191) q[0];
sx q[0];
rz(2.7946081) q[0];
x q[1];
rz(0.25881473) q[2];
sx q[2];
rz(-1.0675028) q[2];
sx q[2];
rz(0.81282114) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3841074) q[1];
sx q[1];
rz(-2.2962484) q[1];
sx q[1];
rz(1.2164335) q[1];
rz(-pi) q[2];
x q[2];
rz(0.21680253) q[3];
sx q[3];
rz(-2.4507782) q[3];
sx q[3];
rz(-1.3135214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.60407698) q[2];
sx q[2];
rz(-0.61415577) q[2];
sx q[2];
rz(2.3094731) q[2];
rz(2.4845947) q[3];
sx q[3];
rz(-1.6712345) q[3];
sx q[3];
rz(-2.7895797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81545365) q[0];
sx q[0];
rz(-0.6093381) q[0];
sx q[0];
rz(-2.4916008) q[0];
rz(-0.12282148) q[1];
sx q[1];
rz(-2.717412) q[1];
sx q[1];
rz(2.4688683) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29826001) q[0];
sx q[0];
rz(-1.5290058) q[0];
sx q[0];
rz(3.1298547) q[0];
rz(-pi) q[1];
rz(-1.2956728) q[2];
sx q[2];
rz(-1.2689991) q[2];
sx q[2];
rz(-1.9375305) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.1993701) q[1];
sx q[1];
rz(-1.3749212) q[1];
sx q[1];
rz(-1.0088527) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0905714) q[3];
sx q[3];
rz(-1.5477763) q[3];
sx q[3];
rz(2.1995403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1136721) q[2];
sx q[2];
rz(-1.547926) q[2];
sx q[2];
rz(2.7552354) q[2];
rz(1.9733285) q[3];
sx q[3];
rz(-1.2315653) q[3];
sx q[3];
rz(1.108235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2439709) q[0];
sx q[0];
rz(-1.4742999) q[0];
sx q[0];
rz(-0.058187159) q[0];
rz(2.8367786) q[1];
sx q[1];
rz(-1.1331646) q[1];
sx q[1];
rz(2.286639) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2475654) q[0];
sx q[0];
rz(-1.5024878) q[0];
sx q[0];
rz(1.5724284) q[0];
rz(2.0570175) q[2];
sx q[2];
rz(-2.2006319) q[2];
sx q[2];
rz(1.6597009) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7185229) q[1];
sx q[1];
rz(-2.3221825) q[1];
sx q[1];
rz(1.2889483) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1596572) q[3];
sx q[3];
rz(-2.5469031) q[3];
sx q[3];
rz(2.2583972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.052224) q[2];
sx q[2];
rz(-1.8171909) q[2];
sx q[2];
rz(1.8734056) q[2];
rz(-0.44102272) q[3];
sx q[3];
rz(-0.62097582) q[3];
sx q[3];
rz(-2.2997901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0190401) q[0];
sx q[0];
rz(-2.1903116) q[0];
sx q[0];
rz(0.71504354) q[0];
rz(-0.41636458) q[1];
sx q[1];
rz(-2.6519897) q[1];
sx q[1];
rz(0.29410902) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9884777) q[0];
sx q[0];
rz(-1.5128969) q[0];
sx q[0];
rz(-0.87429509) q[0];
x q[1];
rz(-0.064925504) q[2];
sx q[2];
rz(-1.7749507) q[2];
sx q[2];
rz(-0.72694187) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3336368) q[1];
sx q[1];
rz(-0.2273493) q[1];
sx q[1];
rz(-2.5264986) q[1];
rz(0.77463051) q[3];
sx q[3];
rz(-2.3916158) q[3];
sx q[3];
rz(0.74912723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2405582) q[2];
sx q[2];
rz(-1.3996539) q[2];
sx q[2];
rz(2.412839) q[2];
rz(-0.48480222) q[3];
sx q[3];
rz(-2.1383643) q[3];
sx q[3];
rz(-0.17759594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6840376) q[0];
sx q[0];
rz(-1.3223248) q[0];
sx q[0];
rz(-2.8628602) q[0];
rz(-3.0740671) q[1];
sx q[1];
rz(-2.4261256) q[1];
sx q[1];
rz(2.6356437) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7161432) q[0];
sx q[0];
rz(-1.4977411) q[0];
sx q[0];
rz(-0.55668932) q[0];
x q[1];
rz(1.672411) q[2];
sx q[2];
rz(-1.7676465) q[2];
sx q[2];
rz(-0.54430279) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7655395) q[1];
sx q[1];
rz(-2.3552822) q[1];
sx q[1];
rz(-2.3885352) q[1];
rz(-pi) q[2];
x q[2];
rz(0.89754126) q[3];
sx q[3];
rz(-1.3607303) q[3];
sx q[3];
rz(-2.1830934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.94023306) q[2];
sx q[2];
rz(-1.6910005) q[2];
sx q[2];
rz(-0.10227164) q[2];
rz(2.8376132) q[3];
sx q[3];
rz(-2.961048) q[3];
sx q[3];
rz(0.48039082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85222307) q[0];
sx q[0];
rz(-0.82813534) q[0];
sx q[0];
rz(0.68036675) q[0];
rz(0.96087372) q[1];
sx q[1];
rz(-1.5545132) q[1];
sx q[1];
rz(0.56112498) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38491098) q[0];
sx q[0];
rz(-1.9327243) q[0];
sx q[0];
rz(-0.2189526) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8041925) q[2];
sx q[2];
rz(-2.4037558) q[2];
sx q[2];
rz(-2.3349175) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1316072) q[1];
sx q[1];
rz(-1.6254043) q[1];
sx q[1];
rz(-2.8427678) q[1];
rz(-0.64021982) q[3];
sx q[3];
rz(-1.9915765) q[3];
sx q[3];
rz(-0.32395485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9095824) q[2];
sx q[2];
rz(-0.97753111) q[2];
sx q[2];
rz(1.5634465) q[2];
rz(-1.0163418) q[3];
sx q[3];
rz(-1.4278102) q[3];
sx q[3];
rz(-2.0411172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4119754) q[0];
sx q[0];
rz(-2.602674) q[0];
sx q[0];
rz(-0.49402657) q[0];
rz(1.5771075) q[1];
sx q[1];
rz(-2.1421075) q[1];
sx q[1];
rz(0.090242537) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9245042) q[0];
sx q[0];
rz(-0.49744931) q[0];
sx q[0];
rz(2.6352232) q[0];
rz(-pi) q[1];
rz(0.87553067) q[2];
sx q[2];
rz(-1.1644496) q[2];
sx q[2];
rz(-1.6900495) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.75587326) q[1];
sx q[1];
rz(-1.8579646) q[1];
sx q[1];
rz(-1.1509291) q[1];
rz(-0.97848864) q[3];
sx q[3];
rz(-0.33708015) q[3];
sx q[3];
rz(-0.85142985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.51284853) q[2];
sx q[2];
rz(-0.57400846) q[2];
sx q[2];
rz(-0.64344704) q[2];
rz(2.511034) q[3];
sx q[3];
rz(-0.75391155) q[3];
sx q[3];
rz(-0.049230922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.028320463) q[0];
sx q[0];
rz(-1.7009108) q[0];
sx q[0];
rz(1.2233618) q[0];
rz(2.2471097) q[1];
sx q[1];
rz(-2.5136785) q[1];
sx q[1];
rz(1.1462513) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.664331) q[0];
sx q[0];
rz(-2.1210175) q[0];
sx q[0];
rz(2.6653637) q[0];
rz(-pi) q[1];
rz(0.38739631) q[2];
sx q[2];
rz(-1.9700389) q[2];
sx q[2];
rz(2.7170167) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.97630298) q[1];
sx q[1];
rz(-2.7600651) q[1];
sx q[1];
rz(-0.69451992) q[1];
x q[2];
rz(-1.9548863) q[3];
sx q[3];
rz(-0.62058228) q[3];
sx q[3];
rz(0.24031249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2724096) q[2];
sx q[2];
rz(-2.2792008) q[2];
sx q[2];
rz(-1.4746846) q[2];
rz(-0.21371755) q[3];
sx q[3];
rz(-1.7361879) q[3];
sx q[3];
rz(2.2739482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2992582) q[0];
sx q[0];
rz(-2.5304351) q[0];
sx q[0];
rz(-0.70511955) q[0];
rz(2.5662388) q[1];
sx q[1];
rz(-1.7333142) q[1];
sx q[1];
rz(0.56952482) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8815959) q[0];
sx q[0];
rz(-1.4896797) q[0];
sx q[0];
rz(-0.077395721) q[0];
x q[1];
rz(1.2766383) q[2];
sx q[2];
rz(-2.0450971) q[2];
sx q[2];
rz(-2.2612342) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3635253) q[1];
sx q[1];
rz(-2.229564) q[1];
sx q[1];
rz(2.3450065) q[1];
rz(-pi) q[2];
rz(1.873718) q[3];
sx q[3];
rz(-2.0407131) q[3];
sx q[3];
rz(-0.63602122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9465955) q[2];
sx q[2];
rz(-2.1985998) q[2];
sx q[2];
rz(-2.1627964) q[2];
rz(2.7106674) q[3];
sx q[3];
rz(-0.48723358) q[3];
sx q[3];
rz(-1.1270969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8579213) q[0];
sx q[0];
rz(-3.1223174) q[0];
sx q[0];
rz(1.4625782) q[0];
rz(2.1025533) q[1];
sx q[1];
rz(-1.3835399) q[1];
sx q[1];
rz(-1.9336112) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6039654) q[0];
sx q[0];
rz(-2.6485291) q[0];
sx q[0];
rz(2.2490671) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1080148) q[2];
sx q[2];
rz(-1.0462282) q[2];
sx q[2];
rz(-1.3474423) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.12700272) q[1];
sx q[1];
rz(-2.5095438) q[1];
sx q[1];
rz(-1.3227425) q[1];
x q[2];
rz(3.0737526) q[3];
sx q[3];
rz(-1.1761936) q[3];
sx q[3];
rz(-2.2369235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.98564467) q[2];
sx q[2];
rz(-1.0040823) q[2];
sx q[2];
rz(-1.3204302) q[2];
rz(2.7909347) q[3];
sx q[3];
rz(-2.3242293) q[3];
sx q[3];
rz(2.5613274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88631267) q[0];
sx q[0];
rz(-1.9552312) q[0];
sx q[0];
rz(2.2807688) q[0];
rz(-1.6571922) q[1];
sx q[1];
rz(-1.3927554) q[1];
sx q[1];
rz(1.4295255) q[1];
rz(-1.6419353) q[2];
sx q[2];
rz(-1.1924469) q[2];
sx q[2];
rz(3.057657) q[2];
rz(-1.2091533) q[3];
sx q[3];
rz(-2.4222838) q[3];
sx q[3];
rz(-2.1771976) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
