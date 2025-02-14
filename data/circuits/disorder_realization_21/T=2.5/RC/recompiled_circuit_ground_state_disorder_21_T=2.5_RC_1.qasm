OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.8706239) q[0];
sx q[0];
rz(-2.5854817) q[0];
sx q[0];
rz(-2.1882353) q[0];
rz(0.019729992) q[1];
sx q[1];
rz(-2.1058197) q[1];
sx q[1];
rz(-2.1929725) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0625437) q[0];
sx q[0];
rz(-1.9808852) q[0];
sx q[0];
rz(2.5586302) q[0];
rz(0.011587338) q[2];
sx q[2];
rz(-0.23749781) q[2];
sx q[2];
rz(1.604014) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.1114914) q[1];
sx q[1];
rz(-1.1426569) q[1];
sx q[1];
rz(0.74049048) q[1];
rz(-1.0179187) q[3];
sx q[3];
rz(-1.8686668) q[3];
sx q[3];
rz(-1.8436933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3775776) q[2];
sx q[2];
rz(-3.0674051) q[2];
sx q[2];
rz(-0.78262502) q[2];
rz(0.11893663) q[3];
sx q[3];
rz(-2.1075893) q[3];
sx q[3];
rz(0.14757601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19621944) q[0];
sx q[0];
rz(-2.0202899) q[0];
sx q[0];
rz(-0.25783208) q[0];
rz(-3.0505772) q[1];
sx q[1];
rz(-2.0423404) q[1];
sx q[1];
rz(-1.6450504) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.634406) q[0];
sx q[0];
rz(-1.0189462) q[0];
sx q[0];
rz(0.080291434) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.37249724) q[2];
sx q[2];
rz(-1.8939233) q[2];
sx q[2];
rz(1.966764) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1196668) q[1];
sx q[1];
rz(-0.41321427) q[1];
sx q[1];
rz(-1.4959123) q[1];
x q[2];
rz(-1.1608382) q[3];
sx q[3];
rz(-2.1968968) q[3];
sx q[3];
rz(-0.61625851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1220793) q[2];
sx q[2];
rz(-1.8668819) q[2];
sx q[2];
rz(0.40461928) q[2];
rz(2.1520065) q[3];
sx q[3];
rz(-1.3000969) q[3];
sx q[3];
rz(2.5185481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0290381) q[0];
sx q[0];
rz(-2.9974388) q[0];
sx q[0];
rz(-0.83830225) q[0];
rz(2.2932032) q[1];
sx q[1];
rz(-1.219607) q[1];
sx q[1];
rz(-2.1489876) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.557195) q[0];
sx q[0];
rz(-2.2879172) q[0];
sx q[0];
rz(-1.756987) q[0];
rz(-pi) q[1];
rz(-2.8072629) q[2];
sx q[2];
rz(-1.5168166) q[2];
sx q[2];
rz(2.8950429) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6684718) q[1];
sx q[1];
rz(-1.8278012) q[1];
sx q[1];
rz(-0.42196749) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8512673) q[3];
sx q[3];
rz(-1.464352) q[3];
sx q[3];
rz(-2.8255879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.8695716) q[2];
sx q[2];
rz(-1.1676936) q[2];
sx q[2];
rz(-1.1248379) q[2];
rz(-2.520842) q[3];
sx q[3];
rz(-0.97095942) q[3];
sx q[3];
rz(3.0264405) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2426185) q[0];
sx q[0];
rz(-1.9816575) q[0];
sx q[0];
rz(-2.9677891) q[0];
rz(-0.68570343) q[1];
sx q[1];
rz(-1.655429) q[1];
sx q[1];
rz(0.83820835) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.099412709) q[0];
sx q[0];
rz(-1.1118708) q[0];
sx q[0];
rz(2.6290429) q[0];
x q[1];
rz(0.72553749) q[2];
sx q[2];
rz(-0.71070403) q[2];
sx q[2];
rz(0.85884554) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.61889583) q[1];
sx q[1];
rz(-2.451183) q[1];
sx q[1];
rz(-2.9179179) q[1];
rz(-pi) q[2];
x q[2];
rz(1.533314) q[3];
sx q[3];
rz(-2.0333383) q[3];
sx q[3];
rz(-1.2087865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.54245) q[2];
sx q[2];
rz(-1.4582381) q[2];
sx q[2];
rz(2.7899) q[2];
rz(-1.9463978) q[3];
sx q[3];
rz(-1.1870793) q[3];
sx q[3];
rz(3.1174507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46566063) q[0];
sx q[0];
rz(-2.6213578) q[0];
sx q[0];
rz(-2.7886673) q[0];
rz(2.832761) q[1];
sx q[1];
rz(-1.0363204) q[1];
sx q[1];
rz(0.028506361) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7362979) q[0];
sx q[0];
rz(-0.60966821) q[0];
sx q[0];
rz(1.8788615) q[0];
rz(-2.7598031) q[2];
sx q[2];
rz(-1.2068527) q[2];
sx q[2];
rz(0.84944968) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5529768) q[1];
sx q[1];
rz(-1.3541095) q[1];
sx q[1];
rz(2.3439247) q[1];
x q[2];
rz(-1.1358374) q[3];
sx q[3];
rz(-1.3490145) q[3];
sx q[3];
rz(-1.3835761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.21294022) q[2];
sx q[2];
rz(-0.063455909) q[2];
sx q[2];
rz(2.1707936) q[2];
rz(1.0468696) q[3];
sx q[3];
rz(-1.4528843) q[3];
sx q[3];
rz(-2.651732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30188072) q[0];
sx q[0];
rz(-1.3588926) q[0];
sx q[0];
rz(-0.090959892) q[0];
rz(1.2184527) q[1];
sx q[1];
rz(-2.8107042) q[1];
sx q[1];
rz(-2.5023696) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6240182) q[0];
sx q[0];
rz(-1.7057452) q[0];
sx q[0];
rz(-2.8132416) q[0];
x q[1];
rz(-0.12219723) q[2];
sx q[2];
rz(-1.4336042) q[2];
sx q[2];
rz(-0.73681632) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6216633) q[1];
sx q[1];
rz(-0.93931336) q[1];
sx q[1];
rz(-2.2261621) q[1];
x q[2];
rz(-0.0019885824) q[3];
sx q[3];
rz(-1.8973777) q[3];
sx q[3];
rz(1.1359362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0794534) q[2];
sx q[2];
rz(-2.1663351) q[2];
sx q[2];
rz(-2.6677168) q[2];
rz(1.3300995) q[3];
sx q[3];
rz(-2.2606943) q[3];
sx q[3];
rz(0.37548319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5652931) q[0];
sx q[0];
rz(-2.0602891) q[0];
sx q[0];
rz(-2.9190049) q[0];
rz(1.0635771) q[1];
sx q[1];
rz(-1.5779326) q[1];
sx q[1];
rz(1.9482013) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.397307) q[0];
sx q[0];
rz(-2.3145876) q[0];
sx q[0];
rz(1.4124407) q[0];
rz(-pi) q[1];
rz(1.3435828) q[2];
sx q[2];
rz(-1.4526723) q[2];
sx q[2];
rz(-0.83641499) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4970253) q[1];
sx q[1];
rz(-2.679832) q[1];
sx q[1];
rz(-2.8341932) q[1];
rz(-pi) q[2];
rz(-1.0289331) q[3];
sx q[3];
rz(-0.30908424) q[3];
sx q[3];
rz(2.2268471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0003164) q[2];
sx q[2];
rz(-0.87339425) q[2];
sx q[2];
rz(-0.6401965) q[2];
rz(0.021942465) q[3];
sx q[3];
rz(-1.4098189) q[3];
sx q[3];
rz(-2.791361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6496277) q[0];
sx q[0];
rz(-1.7438629) q[0];
sx q[0];
rz(-2.6212027) q[0];
rz(0.87977663) q[1];
sx q[1];
rz(-1.2326515) q[1];
sx q[1];
rz(-0.89281503) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2498024) q[0];
sx q[0];
rz(-3.043104) q[0];
sx q[0];
rz(-0.47827025) q[0];
rz(-2.3273349) q[2];
sx q[2];
rz(-0.99254464) q[2];
sx q[2];
rz(1.949343) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8930606) q[1];
sx q[1];
rz(-1.8527781) q[1];
sx q[1];
rz(-2.3470013) q[1];
rz(-pi) q[2];
rz(0.67075394) q[3];
sx q[3];
rz(-0.62707147) q[3];
sx q[3];
rz(2.1334836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2616547) q[2];
sx q[2];
rz(-1.1038021) q[2];
sx q[2];
rz(-3.0547764) q[2];
rz(-0.9451198) q[3];
sx q[3];
rz(-2.7961531) q[3];
sx q[3];
rz(-2.0115578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18411186) q[0];
sx q[0];
rz(-3.0511973) q[0];
sx q[0];
rz(-0.064067319) q[0];
rz(-2.0059026) q[1];
sx q[1];
rz(-1.7013197) q[1];
sx q[1];
rz(0.51220977) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50822608) q[0];
sx q[0];
rz(-2.6056943) q[0];
sx q[0];
rz(0.6750489) q[0];
rz(-1.1697328) q[2];
sx q[2];
rz(-2.3980015) q[2];
sx q[2];
rz(2.0498073) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.87742245) q[1];
sx q[1];
rz(-1.8296731) q[1];
sx q[1];
rz(-1.4822465) q[1];
rz(-pi) q[2];
rz(1.5867611) q[3];
sx q[3];
rz(-1.9796238) q[3];
sx q[3];
rz(0.4086993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3247437) q[2];
sx q[2];
rz(-2.4090359) q[2];
sx q[2];
rz(-2.0762439) q[2];
rz(-1.4893074) q[3];
sx q[3];
rz(-2.1842712) q[3];
sx q[3];
rz(3.0491507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5866933) q[0];
sx q[0];
rz(-0.41189343) q[0];
sx q[0];
rz(-2.8107585) q[0];
rz(1.2384442) q[1];
sx q[1];
rz(-0.60791433) q[1];
sx q[1];
rz(-2.8668561) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82652826) q[0];
sx q[0];
rz(-1.2312268) q[0];
sx q[0];
rz(1.1134321) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0255453) q[2];
sx q[2];
rz(-0.2751285) q[2];
sx q[2];
rz(-2.3809718) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.52636056) q[1];
sx q[1];
rz(-0.52605275) q[1];
sx q[1];
rz(1.8520717) q[1];
x q[2];
rz(-1.6791183) q[3];
sx q[3];
rz(-0.59039718) q[3];
sx q[3];
rz(-0.27449671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9773418) q[2];
sx q[2];
rz(-2.3873316) q[2];
sx q[2];
rz(1.3573307) q[2];
rz(2.6409798) q[3];
sx q[3];
rz(-1.5784135) q[3];
sx q[3];
rz(1.4969426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(1.5464583) q[0];
sx q[0];
rz(-1.5820137) q[0];
sx q[0];
rz(-0.59376846) q[0];
rz(0.068269923) q[1];
sx q[1];
rz(-1.9207813) q[1];
sx q[1];
rz(-3.1178738) q[1];
rz(2.7958128) q[2];
sx q[2];
rz(-1.7263392) q[2];
sx q[2];
rz(1.6172258) q[2];
rz(-2.3292062) q[3];
sx q[3];
rz(-0.25573041) q[3];
sx q[3];
rz(-2.0144016) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
