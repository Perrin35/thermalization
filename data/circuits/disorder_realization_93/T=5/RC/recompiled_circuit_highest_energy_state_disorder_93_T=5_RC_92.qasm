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
rz(2.6858202) q[0];
sx q[0];
rz(-2.0210285) q[0];
sx q[0];
rz(-1.4135345) q[0];
rz(-2.7780374) q[1];
sx q[1];
rz(-1.8283365) q[1];
sx q[1];
rz(-2.8503836) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7257627) q[0];
sx q[0];
rz(-0.32074577) q[0];
sx q[0];
rz(1.6345535) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0179196) q[2];
sx q[2];
rz(-1.9078662) q[2];
sx q[2];
rz(0.45749422) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0800228) q[1];
sx q[1];
rz(-2.1407619) q[1];
sx q[1];
rz(1.7504196) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.81756567) q[3];
sx q[3];
rz(-1.0888087) q[3];
sx q[3];
rz(-2.7369902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.0061965813) q[2];
sx q[2];
rz(-1.5634894) q[2];
sx q[2];
rz(-0.023999365) q[2];
rz(-0.35605797) q[3];
sx q[3];
rz(-0.67432109) q[3];
sx q[3];
rz(-0.02478987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0559219) q[0];
sx q[0];
rz(-2.4503158) q[0];
sx q[0];
rz(-0.63496494) q[0];
rz(0.7629281) q[1];
sx q[1];
rz(-1.678391) q[1];
sx q[1];
rz(-2.2244804) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13697727) q[0];
sx q[0];
rz(-3.1023623) q[0];
sx q[0];
rz(1.0835539) q[0];
x q[1];
rz(-1.00441) q[2];
sx q[2];
rz(-0.76493636) q[2];
sx q[2];
rz(-1.1980008) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.13027993) q[1];
sx q[1];
rz(-0.44751274) q[1];
sx q[1];
rz(1.391973) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.21288721) q[3];
sx q[3];
rz(-0.96159426) q[3];
sx q[3];
rz(1.4573801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1941173) q[2];
sx q[2];
rz(-0.52299356) q[2];
sx q[2];
rz(0.73327649) q[2];
rz(1.0670079) q[3];
sx q[3];
rz(-1.7206444) q[3];
sx q[3];
rz(-0.1799306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5697524) q[0];
sx q[0];
rz(-1.1886007) q[0];
sx q[0];
rz(0.52823129) q[0];
rz(-2.275548) q[1];
sx q[1];
rz(-1.8128315) q[1];
sx q[1];
rz(-2.6285062) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24241867) q[0];
sx q[0];
rz(-1.1326101) q[0];
sx q[0];
rz(0.04695453) q[0];
rz(-pi) q[1];
rz(-1.1703067) q[2];
sx q[2];
rz(-2.4520564) q[2];
sx q[2];
rz(-0.51338085) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9218623) q[1];
sx q[1];
rz(-0.37223909) q[1];
sx q[1];
rz(-1.2138741) q[1];
rz(-pi) q[2];
rz(-1.5422899) q[3];
sx q[3];
rz(-1.0712726) q[3];
sx q[3];
rz(-2.0716425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.036006007) q[2];
sx q[2];
rz(-2.0522108) q[2];
sx q[2];
rz(0.038724381) q[2];
rz(0.45768467) q[3];
sx q[3];
rz(-0.80515146) q[3];
sx q[3];
rz(1.3409486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7078581) q[0];
sx q[0];
rz(-2.7271294) q[0];
sx q[0];
rz(-1.1210972) q[0];
rz(2.3838249) q[1];
sx q[1];
rz(-1.517375) q[1];
sx q[1];
rz(2.4574492) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6586886) q[0];
sx q[0];
rz(-1.7724121) q[0];
sx q[0];
rz(0.828142) q[0];
rz(-pi) q[1];
rz(0.010104601) q[2];
sx q[2];
rz(-1.5779543) q[2];
sx q[2];
rz(-2.7814076) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0348548) q[1];
sx q[1];
rz(-1.8632006) q[1];
sx q[1];
rz(0.3180983) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2627801) q[3];
sx q[3];
rz(-2.4185816) q[3];
sx q[3];
rz(0.77817749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6571558) q[2];
sx q[2];
rz(-1.6146722) q[2];
sx q[2];
rz(0.26051513) q[2];
rz(-0.83812964) q[3];
sx q[3];
rz(-1.9939491) q[3];
sx q[3];
rz(-3.1269791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92358661) q[0];
sx q[0];
rz(-1.3777233) q[0];
sx q[0];
rz(2.6309784) q[0];
rz(-1.249373) q[1];
sx q[1];
rz(-1.1621954) q[1];
sx q[1];
rz(-1.6872663) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73390612) q[0];
sx q[0];
rz(-0.21557325) q[0];
sx q[0];
rz(-3.0093497) q[0];
x q[1];
rz(0.32938214) q[2];
sx q[2];
rz(-2.9940146) q[2];
sx q[2];
rz(-1.3775228) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.37794681) q[1];
sx q[1];
rz(-2.1140522) q[1];
sx q[1];
rz(1.2681947) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9972059) q[3];
sx q[3];
rz(-1.4098415) q[3];
sx q[3];
rz(-0.99100366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.55383596) q[2];
sx q[2];
rz(-0.38148701) q[2];
sx q[2];
rz(-0.39056632) q[2];
rz(-2.8824814) q[3];
sx q[3];
rz(-1.6389537) q[3];
sx q[3];
rz(0.34671569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7873586) q[0];
sx q[0];
rz(-0.43208313) q[0];
sx q[0];
rz(2.4969192) q[0];
rz(-1.8395754) q[1];
sx q[1];
rz(-1.6141067) q[1];
sx q[1];
rz(-0.34861809) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7946588) q[0];
sx q[0];
rz(-2.0290792) q[0];
sx q[0];
rz(-2.6285104) q[0];
rz(-pi) q[1];
rz(-2.5041144) q[2];
sx q[2];
rz(-0.63796746) q[2];
sx q[2];
rz(1.5359985) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8585414) q[1];
sx q[1];
rz(-1.9894454) q[1];
sx q[1];
rz(0.323723) q[1];
x q[2];
rz(-3.0246234) q[3];
sx q[3];
rz(-1.5364963) q[3];
sx q[3];
rz(-0.12730612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1031441) q[2];
sx q[2];
rz(-0.95547533) q[2];
sx q[2];
rz(-3.1032584) q[2];
rz(-2.1740055) q[3];
sx q[3];
rz(-1.4097694) q[3];
sx q[3];
rz(-2.7819395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2751145) q[0];
sx q[0];
rz(-2.622128) q[0];
sx q[0];
rz(3.0679829) q[0];
rz(1.4069125) q[1];
sx q[1];
rz(-2.584447) q[1];
sx q[1];
rz(0.73307347) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.276537) q[0];
sx q[0];
rz(-0.77651327) q[0];
sx q[0];
rz(-0.76188253) q[0];
x q[1];
rz(-2.5400852) q[2];
sx q[2];
rz(-0.93346338) q[2];
sx q[2];
rz(1.9236227) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.16790529) q[1];
sx q[1];
rz(-1.2406772) q[1];
sx q[1];
rz(0.37574212) q[1];
rz(-0.81871512) q[3];
sx q[3];
rz(-0.33740852) q[3];
sx q[3];
rz(-0.81160566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1072032) q[2];
sx q[2];
rz(-1.723039) q[2];
sx q[2];
rz(-3.1305195) q[2];
rz(0.30558807) q[3];
sx q[3];
rz(-1.1253076) q[3];
sx q[3];
rz(-1.8886458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8082064) q[0];
sx q[0];
rz(-0.60788637) q[0];
sx q[0];
rz(0.73390865) q[0];
rz(-0.58087307) q[1];
sx q[1];
rz(-2.1323233) q[1];
sx q[1];
rz(-3.0427921) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2487468) q[0];
sx q[0];
rz(-2.1552784) q[0];
sx q[0];
rz(-2.2131323) q[0];
rz(-pi) q[1];
x q[1];
rz(0.67811857) q[2];
sx q[2];
rz(-2.0653915) q[2];
sx q[2];
rz(-1.8910318) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.1360201) q[1];
sx q[1];
rz(-0.51678777) q[1];
sx q[1];
rz(2.6173475) q[1];
rz(-pi) q[2];
rz(1.6510165) q[3];
sx q[3];
rz(-1.7003683) q[3];
sx q[3];
rz(-3.0440273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.29844478) q[2];
sx q[2];
rz(-1.8567905) q[2];
sx q[2];
rz(-1.0996381) q[2];
rz(-3.1067276) q[3];
sx q[3];
rz(-2.3979082) q[3];
sx q[3];
rz(0.59665027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.027503969) q[0];
sx q[0];
rz(-1.9845668) q[0];
sx q[0];
rz(-1.9644568) q[0];
rz(-1.6674532) q[1];
sx q[1];
rz(-0.93282229) q[1];
sx q[1];
rz(1.4449545) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6072877) q[0];
sx q[0];
rz(-1.289708) q[0];
sx q[0];
rz(2.7318153) q[0];
rz(0.99924008) q[2];
sx q[2];
rz(-2.2543467) q[2];
sx q[2];
rz(-0.37887606) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7921389) q[1];
sx q[1];
rz(-1.0568403) q[1];
sx q[1];
rz(0.07899125) q[1];
x q[2];
rz(-2.2463465) q[3];
sx q[3];
rz(-1.352868) q[3];
sx q[3];
rz(-2.982614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.16704796) q[2];
sx q[2];
rz(-1.5571152) q[2];
sx q[2];
rz(2.9696999) q[2];
rz(0.78957549) q[3];
sx q[3];
rz(-2.004345) q[3];
sx q[3];
rz(-1.8062228) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0563141) q[0];
sx q[0];
rz(-0.21509898) q[0];
sx q[0];
rz(-2.6000182) q[0];
rz(1.1594634) q[1];
sx q[1];
rz(-1.3075202) q[1];
sx q[1];
rz(0.83311876) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.706382) q[0];
sx q[0];
rz(-0.66775371) q[0];
sx q[0];
rz(0.45489648) q[0];
rz(-1.2255185) q[2];
sx q[2];
rz(-2.08689) q[2];
sx q[2];
rz(-0.47613019) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7685582) q[1];
sx q[1];
rz(-1.3443861) q[1];
sx q[1];
rz(-0.86707244) q[1];
x q[2];
rz(1.5686137) q[3];
sx q[3];
rz(-2.4139521) q[3];
sx q[3];
rz(-1.1436056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.1200166) q[2];
sx q[2];
rz(-0.77777255) q[2];
sx q[2];
rz(-0.94338256) q[2];
rz(-2.0639482) q[3];
sx q[3];
rz(-1.5871983) q[3];
sx q[3];
rz(0.32850346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
rz(-0.42644603) q[0];
sx q[0];
rz(-1.0363415) q[0];
sx q[0];
rz(2.8847726) q[0];
rz(2.9019451) q[1];
sx q[1];
rz(-0.9050723) q[1];
sx q[1];
rz(-1.1368652) q[1];
rz(0.73020936) q[2];
sx q[2];
rz(-2.3013023) q[2];
sx q[2];
rz(-0.58462044) q[2];
rz(-1.3171468) q[3];
sx q[3];
rz(-1.2696878) q[3];
sx q[3];
rz(2.9326143) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
