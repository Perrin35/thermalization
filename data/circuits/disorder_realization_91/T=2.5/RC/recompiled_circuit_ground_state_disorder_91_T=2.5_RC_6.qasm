OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.33136) q[0];
sx q[0];
rz(-1.2547837) q[0];
sx q[0];
rz(-2.4956508) q[0];
rz(1.6826001) q[1];
sx q[1];
rz(-3.0134633) q[1];
sx q[1];
rz(0.66817862) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2006045) q[0];
sx q[0];
rz(-1.4580237) q[0];
sx q[0];
rz(2.8975945) q[0];
rz(-pi) q[1];
rz(0.3273079) q[2];
sx q[2];
rz(-2.1919554) q[2];
sx q[2];
rz(-2.6747764) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9757802) q[1];
sx q[1];
rz(-2.4620612) q[1];
sx q[1];
rz(0.91109101) q[1];
rz(-0.48623881) q[3];
sx q[3];
rz(-0.91043962) q[3];
sx q[3];
rz(2.5442991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.51840034) q[2];
sx q[2];
rz(-0.65541583) q[2];
sx q[2];
rz(-2.0584959) q[2];
rz(-0.25534233) q[3];
sx q[3];
rz(-0.89699236) q[3];
sx q[3];
rz(1.159509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1622247) q[0];
sx q[0];
rz(-0.14142445) q[0];
sx q[0];
rz(-3.1080143) q[0];
rz(-2.0032739) q[1];
sx q[1];
rz(-0.99457026) q[1];
sx q[1];
rz(2.9046362) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19582307) q[0];
sx q[0];
rz(-1.3422215) q[0];
sx q[0];
rz(-1.7128471) q[0];
rz(2.9541624) q[2];
sx q[2];
rz(-2.0399562) q[2];
sx q[2];
rz(-0.7721061) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0396871) q[1];
sx q[1];
rz(-2.9475005) q[1];
sx q[1];
rz(-1.4386574) q[1];
rz(-0.44741607) q[3];
sx q[3];
rz(-1.0378812) q[3];
sx q[3];
rz(0.79870236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.61243764) q[2];
sx q[2];
rz(-1.0470904) q[2];
sx q[2];
rz(3.0212413) q[2];
rz(0.58593166) q[3];
sx q[3];
rz(-1.3804932) q[3];
sx q[3];
rz(1.8732635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2633857) q[0];
sx q[0];
rz(-0.96579856) q[0];
sx q[0];
rz(2.1534488) q[0];
rz(-1.6506317) q[1];
sx q[1];
rz(-2.4246876) q[1];
sx q[1];
rz(1.5879226) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7415888) q[0];
sx q[0];
rz(-0.015946139) q[0];
sx q[0];
rz(3.0724597) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0071909) q[2];
sx q[2];
rz(-1.5024868) q[2];
sx q[2];
rz(1.3775064) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2026489) q[1];
sx q[1];
rz(-1.5420785) q[1];
sx q[1];
rz(1.3095301) q[1];
rz(2.5293143) q[3];
sx q[3];
rz(-0.66198549) q[3];
sx q[3];
rz(0.16890165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.83739088) q[2];
sx q[2];
rz(-1.7981671) q[2];
sx q[2];
rz(3.0710132) q[2];
rz(0.48921674) q[3];
sx q[3];
rz(-2.1507806) q[3];
sx q[3];
rz(-2.9941471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9551142) q[0];
sx q[0];
rz(-0.69545737) q[0];
sx q[0];
rz(1.7406933) q[0];
rz(-0.1700302) q[1];
sx q[1];
rz(-1.731512) q[1];
sx q[1];
rz(1.2331351) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2182541) q[0];
sx q[0];
rz(-1.4460576) q[0];
sx q[0];
rz(1.1472923) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.71738404) q[2];
sx q[2];
rz(-2.0307433) q[2];
sx q[2];
rz(-1.1954952) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4608232) q[1];
sx q[1];
rz(-1.5857547) q[1];
sx q[1];
rz(-1.8504924) q[1];
rz(-0.024007576) q[3];
sx q[3];
rz(-2.6769961) q[3];
sx q[3];
rz(-0.28991227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.79973811) q[2];
sx q[2];
rz(-1.7065455) q[2];
sx q[2];
rz(-0.62961659) q[2];
rz(-0.13684212) q[3];
sx q[3];
rz(-0.62961951) q[3];
sx q[3];
rz(1.4153882) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8585696) q[0];
sx q[0];
rz(-0.50823277) q[0];
sx q[0];
rz(-3.0392905) q[0];
rz(-1.6434068) q[1];
sx q[1];
rz(-1.3003277) q[1];
sx q[1];
rz(2.7844875) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67960167) q[0];
sx q[0];
rz(-1.7850189) q[0];
sx q[0];
rz(-0.22837436) q[0];
x q[1];
rz(-0.12320766) q[2];
sx q[2];
rz(-1.5438617) q[2];
sx q[2];
rz(-1.3696485) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9041464) q[1];
sx q[1];
rz(-1.0654628) q[1];
sx q[1];
rz(1.2282441) q[1];
x q[2];
rz(0.83729845) q[3];
sx q[3];
rz(-2.5243763) q[3];
sx q[3];
rz(1.069998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6685278) q[2];
sx q[2];
rz(-1.4184971) q[2];
sx q[2];
rz(1.3225887) q[2];
rz(-1.708301) q[3];
sx q[3];
rz(-1.3168443) q[3];
sx q[3];
rz(-2.9776261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.598031) q[0];
sx q[0];
rz(-0.99995166) q[0];
sx q[0];
rz(-1.438197) q[0];
rz(1.2403129) q[1];
sx q[1];
rz(-2.4867609) q[1];
sx q[1];
rz(-1.7887438) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.529218) q[0];
sx q[0];
rz(-1.9758401) q[0];
sx q[0];
rz(-0.10251001) q[0];
rz(-0.5261658) q[2];
sx q[2];
rz(-1.4246829) q[2];
sx q[2];
rz(2.7046159) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2861915) q[1];
sx q[1];
rz(-1.0637861) q[1];
sx q[1];
rz(0.42393522) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6840965) q[3];
sx q[3];
rz(-1.9637412) q[3];
sx q[3];
rz(-1.9696152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9271586) q[2];
sx q[2];
rz(-2.7775601) q[2];
sx q[2];
rz(-0.39043179) q[2];
rz(1.8292142) q[3];
sx q[3];
rz(-2.3305011) q[3];
sx q[3];
rz(-0.81982476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.696233) q[0];
sx q[0];
rz(-0.92388988) q[0];
sx q[0];
rz(1.4129289) q[0];
rz(-1.6784809) q[1];
sx q[1];
rz(-1.585958) q[1];
sx q[1];
rz(2.5020592) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93776751) q[0];
sx q[0];
rz(-0.97903189) q[0];
sx q[0];
rz(-0.45543619) q[0];
rz(0.15070559) q[2];
sx q[2];
rz(-1.6455212) q[2];
sx q[2];
rz(-3.0091803) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.29296103) q[1];
sx q[1];
rz(-0.94215122) q[1];
sx q[1];
rz(1.2413003) q[1];
rz(-pi) q[2];
rz(-1.8516435) q[3];
sx q[3];
rz(-2.7327545) q[3];
sx q[3];
rz(-2.051681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9245236) q[2];
sx q[2];
rz(-1.8127433) q[2];
sx q[2];
rz(2.6893943) q[2];
rz(0.2229812) q[3];
sx q[3];
rz(-2.860234) q[3];
sx q[3];
rz(0.15534672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13101354) q[0];
sx q[0];
rz(-2.2192945) q[0];
sx q[0];
rz(0.16192326) q[0];
rz(0.34795347) q[1];
sx q[1];
rz(-1.8674928) q[1];
sx q[1];
rz(0.23439342) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97353712) q[0];
sx q[0];
rz(-2.148365) q[0];
sx q[0];
rz(0.67649354) q[0];
rz(-pi) q[1];
rz(-2.4215464) q[2];
sx q[2];
rz(-1.3338425) q[2];
sx q[2];
rz(0.64575486) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.80658) q[1];
sx q[1];
rz(-1.7579105) q[1];
sx q[1];
rz(0.92430618) q[1];
rz(-pi) q[2];
rz(0.23080821) q[3];
sx q[3];
rz(-2.857702) q[3];
sx q[3];
rz(2.6687255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8552385) q[2];
sx q[2];
rz(-2.0704634) q[2];
sx q[2];
rz(-2.5174649) q[2];
rz(-0.74328077) q[3];
sx q[3];
rz(-2.1478839) q[3];
sx q[3];
rz(0.68964094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8648935) q[0];
sx q[0];
rz(-1.445329) q[0];
sx q[0];
rz(-2.0546761) q[0];
rz(1.2326321) q[1];
sx q[1];
rz(-1.107736) q[1];
sx q[1];
rz(1.8392275) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83643276) q[0];
sx q[0];
rz(-1.4576685) q[0];
sx q[0];
rz(2.9311084) q[0];
rz(0.97424284) q[2];
sx q[2];
rz(-1.8678209) q[2];
sx q[2];
rz(3.0926306) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2375396) q[1];
sx q[1];
rz(-1.9642648) q[1];
sx q[1];
rz(0.57668011) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5518673) q[3];
sx q[3];
rz(-0.94784465) q[3];
sx q[3];
rz(-1.2894693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.46633259) q[2];
sx q[2];
rz(-1.8856498) q[2];
sx q[2];
rz(-0.41435286) q[2];
rz(-0.77053344) q[3];
sx q[3];
rz(-1.7194175) q[3];
sx q[3];
rz(-2.2079302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3122124) q[0];
sx q[0];
rz(-2.3214564) q[0];
sx q[0];
rz(2.6729551) q[0];
rz(-1.8679484) q[1];
sx q[1];
rz(-1.2761152) q[1];
sx q[1];
rz(-2.830107) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0429223) q[0];
sx q[0];
rz(-2.2880974) q[0];
sx q[0];
rz(0.076626549) q[0];
rz(1.6577254) q[2];
sx q[2];
rz(-1.6875899) q[2];
sx q[2];
rz(0.96125666) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2356687) q[1];
sx q[1];
rz(-2.0579268) q[1];
sx q[1];
rz(0.052799932) q[1];
rz(-0.52759513) q[3];
sx q[3];
rz(-2.0012534) q[3];
sx q[3];
rz(2.2077219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1062539) q[2];
sx q[2];
rz(-1.1941348) q[2];
sx q[2];
rz(1.3930456) q[2];
rz(-2.2140908) q[3];
sx q[3];
rz(-1.564097) q[3];
sx q[3];
rz(-0.86021304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8504234) q[0];
sx q[0];
rz(-2.8130154) q[0];
sx q[0];
rz(-1.6541506) q[0];
rz(-0.18118478) q[1];
sx q[1];
rz(-2.1962427) q[1];
sx q[1];
rz(-0.93999351) q[1];
rz(1.0370902) q[2];
sx q[2];
rz(-1.4363534) q[2];
sx q[2];
rz(-1.7912122) q[2];
rz(-0.84862205) q[3];
sx q[3];
rz(-0.80905882) q[3];
sx q[3];
rz(-0.73425135) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
