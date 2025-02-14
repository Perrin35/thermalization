OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.085091703) q[0];
sx q[0];
rz(-0.29274517) q[0];
sx q[0];
rz(-0.91896397) q[0];
rz(2.7594944) q[1];
sx q[1];
rz(-0.16799071) q[1];
sx q[1];
rz(-1.1408495) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1534682) q[0];
sx q[0];
rz(-1.690295) q[0];
sx q[0];
rz(1.7480127) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8542641) q[2];
sx q[2];
rz(-1.8745443) q[2];
sx q[2];
rz(1.2992573) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.1124494) q[1];
sx q[1];
rz(-1.1894798) q[1];
sx q[1];
rz(1.7138395) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.010072) q[3];
sx q[3];
rz(-1.4460996) q[3];
sx q[3];
rz(1.6623868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4484278) q[2];
sx q[2];
rz(-1.0591256) q[2];
sx q[2];
rz(-1.9654174) q[2];
rz(-3.0991992) q[3];
sx q[3];
rz(-1.3375125) q[3];
sx q[3];
rz(0.5347518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46928826) q[0];
sx q[0];
rz(-0.35457087) q[0];
sx q[0];
rz(0.18445036) q[0];
rz(-0.09672673) q[1];
sx q[1];
rz(-1.5238785) q[1];
sx q[1];
rz(1.9116037) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7589129) q[0];
sx q[0];
rz(-1.935674) q[0];
sx q[0];
rz(1.7880102) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0646992) q[2];
sx q[2];
rz(-2.580759) q[2];
sx q[2];
rz(-2.4397813) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9877731) q[1];
sx q[1];
rz(-1.7708988) q[1];
sx q[1];
rz(-2.8994183) q[1];
rz(0.038944728) q[3];
sx q[3];
rz(-1.9488261) q[3];
sx q[3];
rz(1.8824487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7868598) q[2];
sx q[2];
rz(-2.730098) q[2];
sx q[2];
rz(3.1301609) q[2];
rz(1.3139906) q[3];
sx q[3];
rz(-1.7766137) q[3];
sx q[3];
rz(0.7029117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74362022) q[0];
sx q[0];
rz(-3.008606) q[0];
sx q[0];
rz(-2.5308894) q[0];
rz(-0.01344219) q[1];
sx q[1];
rz(-1.8553269) q[1];
sx q[1];
rz(-0.59649831) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4134408) q[0];
sx q[0];
rz(-2.1755009) q[0];
sx q[0];
rz(0.20264969) q[0];
rz(-2.1488068) q[2];
sx q[2];
rz(-1.8278549) q[2];
sx q[2];
rz(-1.3959194) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.96287476) q[1];
sx q[1];
rz(-1.2591919) q[1];
sx q[1];
rz(0.7473263) q[1];
rz(-pi) q[2];
rz(3.1085988) q[3];
sx q[3];
rz(-1.0171431) q[3];
sx q[3];
rz(-0.43137303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.67768031) q[2];
sx q[2];
rz(-1.8296506) q[2];
sx q[2];
rz(-0.00027351969) q[2];
rz(1.6655946) q[3];
sx q[3];
rz(-2.4725584) q[3];
sx q[3];
rz(-0.10703787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45604712) q[0];
sx q[0];
rz(-0.16655971) q[0];
sx q[0];
rz(1.1628994) q[0];
rz(-2.1641425) q[1];
sx q[1];
rz(-1.6012499) q[1];
sx q[1];
rz(1.5553442) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88240963) q[0];
sx q[0];
rz(-1.5804187) q[0];
sx q[0];
rz(-3.1338047) q[0];
rz(-pi) q[1];
rz(-2.3562217) q[2];
sx q[2];
rz(-1.5391239) q[2];
sx q[2];
rz(-0.23508628) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0229637) q[1];
sx q[1];
rz(-1.5563772) q[1];
sx q[1];
rz(-0.49683907) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0647207) q[3];
sx q[3];
rz(-1.7887497) q[3];
sx q[3];
rz(-2.2214784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1996475) q[2];
sx q[2];
rz(-0.5117828) q[2];
sx q[2];
rz(-2.9072348) q[2];
rz(-2.6394081) q[3];
sx q[3];
rz(-2.1000704) q[3];
sx q[3];
rz(-2.485399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.31768826) q[0];
sx q[0];
rz(-2.3264139) q[0];
sx q[0];
rz(-1.3932047) q[0];
rz(2.4488917) q[1];
sx q[1];
rz(-1.0857948) q[1];
sx q[1];
rz(2.0511973) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22995725) q[0];
sx q[0];
rz(-1.3595694) q[0];
sx q[0];
rz(-2.1043572) q[0];
rz(2.1887921) q[2];
sx q[2];
rz(-2.0119466) q[2];
sx q[2];
rz(1.49681) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.322256) q[1];
sx q[1];
rz(-2.6988479) q[1];
sx q[1];
rz(-1.2135452) q[1];
rz(-2.9850335) q[3];
sx q[3];
rz(-1.0362175) q[3];
sx q[3];
rz(-1.8434356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.56110567) q[2];
sx q[2];
rz(-1.1539536) q[2];
sx q[2];
rz(-2.6055276) q[2];
rz(-0.29340336) q[3];
sx q[3];
rz(-1.2111726) q[3];
sx q[3];
rz(-2.6611879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2368161) q[0];
sx q[0];
rz(-0.95229709) q[0];
sx q[0];
rz(0.099040898) q[0];
rz(1.0320484) q[1];
sx q[1];
rz(-0.85056225) q[1];
sx q[1];
rz(1.7031857) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28171646) q[0];
sx q[0];
rz(-1.8252354) q[0];
sx q[0];
rz(-0.86812015) q[0];
rz(-pi) q[1];
rz(-0.2834592) q[2];
sx q[2];
rz(-1.4999031) q[2];
sx q[2];
rz(2.9170582) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.21675303) q[1];
sx q[1];
rz(-0.79063639) q[1];
sx q[1];
rz(0.032921493) q[1];
rz(0.29797192) q[3];
sx q[3];
rz(-1.4557585) q[3];
sx q[3];
rz(-1.486766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1556039) q[2];
sx q[2];
rz(-2.6376548) q[2];
sx q[2];
rz(-2.3206553) q[2];
rz(-1.692903) q[3];
sx q[3];
rz(-2.4235348) q[3];
sx q[3];
rz(1.4887571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89707017) q[0];
sx q[0];
rz(-2.2269766) q[0];
sx q[0];
rz(2.6389417) q[0];
rz(-1.9082759) q[1];
sx q[1];
rz(-2.2063467) q[1];
sx q[1];
rz(0.9800235) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9217958) q[0];
sx q[0];
rz(-1.8077407) q[0];
sx q[0];
rz(1.8264997) q[0];
x q[1];
rz(-2.365391) q[2];
sx q[2];
rz(-1.3598249) q[2];
sx q[2];
rz(0.38061505) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7777694) q[1];
sx q[1];
rz(-0.8447434) q[1];
sx q[1];
rz(-1.8822303) q[1];
x q[2];
rz(0.84021826) q[3];
sx q[3];
rz(-1.7257236) q[3];
sx q[3];
rz(-1.2511826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.22214733) q[2];
sx q[2];
rz(-1.8541226) q[2];
sx q[2];
rz(1.3745098) q[2];
rz(1.3162656) q[3];
sx q[3];
rz(-2.2856183) q[3];
sx q[3];
rz(0.98178664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.114349) q[0];
sx q[0];
rz(-0.39334941) q[0];
sx q[0];
rz(2.223176) q[0];
rz(-0.20485993) q[1];
sx q[1];
rz(-2.0826976) q[1];
sx q[1];
rz(1.6580261) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.525104) q[0];
sx q[0];
rz(-1.2813998) q[0];
sx q[0];
rz(-0.15286907) q[0];
rz(-1.4409562) q[2];
sx q[2];
rz(-2.3774638) q[2];
sx q[2];
rz(1.6418599) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.601348) q[1];
sx q[1];
rz(-1.6588604) q[1];
sx q[1];
rz(0.87012) q[1];
rz(-3.0851787) q[3];
sx q[3];
rz(-2.4216363) q[3];
sx q[3];
rz(2.3145793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.001699) q[2];
sx q[2];
rz(-2.6504982) q[2];
sx q[2];
rz(-1.8074544) q[2];
rz(-2.259518) q[3];
sx q[3];
rz(-0.69552723) q[3];
sx q[3];
rz(1.849256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.081721574) q[0];
sx q[0];
rz(-0.43161714) q[0];
sx q[0];
rz(0.01817848) q[0];
rz(0.40953088) q[1];
sx q[1];
rz(-1.7117701) q[1];
sx q[1];
rz(-0.10083625) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1153959) q[0];
sx q[0];
rz(-1.6053204) q[0];
sx q[0];
rz(-2.5071397) q[0];
rz(-pi) q[1];
rz(0.11282044) q[2];
sx q[2];
rz(-2.718528) q[2];
sx q[2];
rz(-3.0213838) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8162569) q[1];
sx q[1];
rz(-2.6995669) q[1];
sx q[1];
rz(-2.1620552) q[1];
x q[2];
rz(-1.3235247) q[3];
sx q[3];
rz(-1.6880638) q[3];
sx q[3];
rz(-0.7922678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.33710256) q[2];
sx q[2];
rz(-0.10919315) q[2];
sx q[2];
rz(1.2479372) q[2];
rz(-1.2248056) q[3];
sx q[3];
rz(-2.2962511) q[3];
sx q[3];
rz(-3.0758408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9058022) q[0];
sx q[0];
rz(-1.6658655) q[0];
sx q[0];
rz(2.8210848) q[0];
rz(0.39101741) q[1];
sx q[1];
rz(-2.0000439) q[1];
sx q[1];
rz(-1.549622) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.029767903) q[0];
sx q[0];
rz(-2.8643908) q[0];
sx q[0];
rz(-0.76407822) q[0];
rz(-pi) q[1];
rz(2.6274849) q[2];
sx q[2];
rz(-1.4288581) q[2];
sx q[2];
rz(2.1588003) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.446881) q[1];
sx q[1];
rz(-1.2031735) q[1];
sx q[1];
rz(-0.093631677) q[1];
x q[2];
rz(-2.7192468) q[3];
sx q[3];
rz(-1.5776488) q[3];
sx q[3];
rz(0.50240483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0512507) q[2];
sx q[2];
rz(-1.0498468) q[2];
sx q[2];
rz(-2.7412097) q[2];
rz(2.9054437) q[3];
sx q[3];
rz(-0.103424) q[3];
sx q[3];
rz(-2.1307438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.878933) q[0];
sx q[0];
rz(-2.6317609) q[0];
sx q[0];
rz(1.2608933) q[0];
rz(-0.48514584) q[1];
sx q[1];
rz(-0.79882516) q[1];
sx q[1];
rz(-0.47732236) q[1];
rz(-1.7791228) q[2];
sx q[2];
rz(-1.4478417) q[2];
sx q[2];
rz(1.7461591) q[2];
rz(0.063592576) q[3];
sx q[3];
rz(-1.7983754) q[3];
sx q[3];
rz(1.6354431) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
