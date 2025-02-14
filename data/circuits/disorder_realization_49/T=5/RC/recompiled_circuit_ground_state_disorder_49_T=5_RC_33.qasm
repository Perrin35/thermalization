OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.9743118) q[0];
sx q[0];
rz(-0.33742961) q[0];
sx q[0];
rz(0.20198241) q[0];
rz(2.1972411) q[1];
sx q[1];
rz(2.386932) q[1];
sx q[1];
rz(7.9805482) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7118397) q[0];
sx q[0];
rz(-2.0637207) q[0];
sx q[0];
rz(-1.809354) q[0];
rz(-3.0336122) q[2];
sx q[2];
rz(-1.3657346) q[2];
sx q[2];
rz(-2.972796) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3814427) q[1];
sx q[1];
rz(-1.1007778) q[1];
sx q[1];
rz(2.2862741) q[1];
x q[2];
rz(-1.6891278) q[3];
sx q[3];
rz(-0.36709309) q[3];
sx q[3];
rz(1.577842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8586388) q[2];
sx q[2];
rz(-0.084246548) q[2];
sx q[2];
rz(-1.5343182) q[2];
rz(-2.9547847) q[3];
sx q[3];
rz(-0.45972937) q[3];
sx q[3];
rz(-1.4927347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9635791) q[0];
sx q[0];
rz(-2.3602965) q[0];
sx q[0];
rz(2.436893) q[0];
rz(-1.0164545) q[1];
sx q[1];
rz(-1.1926788) q[1];
sx q[1];
rz(-1.9889529) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0148852) q[0];
sx q[0];
rz(-2.5683157) q[0];
sx q[0];
rz(-1.8886198) q[0];
rz(-pi) q[1];
rz(0.57781685) q[2];
sx q[2];
rz(-2.0426828) q[2];
sx q[2];
rz(0.058171169) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0699393) q[1];
sx q[1];
rz(-2.3333683) q[1];
sx q[1];
rz(0.83856844) q[1];
rz(-pi) q[2];
rz(-0.14039881) q[3];
sx q[3];
rz(-0.75274668) q[3];
sx q[3];
rz(-1.5206159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.54027259) q[2];
sx q[2];
rz(-0.86869621) q[2];
sx q[2];
rz(-1.3199838) q[2];
rz(3.0705304) q[3];
sx q[3];
rz(-0.080852121) q[3];
sx q[3];
rz(-2.0474056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1048626) q[0];
sx q[0];
rz(-0.34958378) q[0];
sx q[0];
rz(-2.2678251) q[0];
rz(1.8050487) q[1];
sx q[1];
rz(-0.68160325) q[1];
sx q[1];
rz(-0.85493404) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23531518) q[0];
sx q[0];
rz(-1.3585257) q[0];
sx q[0];
rz(-0.92001037) q[0];
rz(-pi) q[1];
rz(0.94742355) q[2];
sx q[2];
rz(-2.6389942) q[2];
sx q[2];
rz(1.1029152) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.47006059) q[1];
sx q[1];
rz(-1.6796659) q[1];
sx q[1];
rz(-2.837869) q[1];
x q[2];
rz(-1.072238) q[3];
sx q[3];
rz(-2.4756458) q[3];
sx q[3];
rz(1.3897105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9285589) q[2];
sx q[2];
rz(-2.010767) q[2];
sx q[2];
rz(-2.0862759) q[2];
rz(-0.94163752) q[3];
sx q[3];
rz(-1.0712653) q[3];
sx q[3];
rz(2.1987703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.4554491) q[0];
sx q[0];
rz(-2.2930155) q[0];
sx q[0];
rz(2.4461179) q[0];
rz(-1.4675568) q[1];
sx q[1];
rz(-1.2662042) q[1];
sx q[1];
rz(3.0470336) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1098925) q[0];
sx q[0];
rz(-1.5227093) q[0];
sx q[0];
rz(-1.2926433) q[0];
x q[1];
rz(0.90717051) q[2];
sx q[2];
rz(-2.2988875) q[2];
sx q[2];
rz(-2.0913194) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.91969216) q[1];
sx q[1];
rz(-1.8522433) q[1];
sx q[1];
rz(-2.9355461) q[1];
x q[2];
rz(-0.88282449) q[3];
sx q[3];
rz(-1.7738288) q[3];
sx q[3];
rz(2.1830851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.62119421) q[2];
sx q[2];
rz(-0.79402557) q[2];
sx q[2];
rz(-2.3801079) q[2];
rz(-2.9705808) q[3];
sx q[3];
rz(-1.2820425) q[3];
sx q[3];
rz(0.096693501) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3068202) q[0];
sx q[0];
rz(-0.67245317) q[0];
sx q[0];
rz(0.49224725) q[0];
rz(-1.3193839) q[1];
sx q[1];
rz(-1.7183813) q[1];
sx q[1];
rz(-0.48670235) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20535417) q[0];
sx q[0];
rz(-1.9689318) q[0];
sx q[0];
rz(1.7004844) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9453508) q[2];
sx q[2];
rz(-0.83714467) q[2];
sx q[2];
rz(-2.9840699) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0472187) q[1];
sx q[1];
rz(-1.1446806) q[1];
sx q[1];
rz(2.7414118) q[1];
x q[2];
rz(-0.23426849) q[3];
sx q[3];
rz(-2.0741012) q[3];
sx q[3];
rz(0.64555321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.95185602) q[2];
sx q[2];
rz(-2.2978013) q[2];
sx q[2];
rz(-2.0437415) q[2];
rz(-2.3073933) q[3];
sx q[3];
rz(-2.4662377) q[3];
sx q[3];
rz(1.437291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-0.51682353) q[0];
sx q[0];
rz(-0.56684816) q[0];
sx q[0];
rz(0.045850642) q[0];
rz(0.74371964) q[1];
sx q[1];
rz(-0.36879483) q[1];
sx q[1];
rz(0.060294453) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3876095) q[0];
sx q[0];
rz(-1.2112677) q[0];
sx q[0];
rz(2.6922168) q[0];
rz(-2.7457775) q[2];
sx q[2];
rz(-2.2278053) q[2];
sx q[2];
rz(0.78570494) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0452779) q[1];
sx q[1];
rz(-1.8966676) q[1];
sx q[1];
rz(2.9067944) q[1];
rz(-0.095796776) q[3];
sx q[3];
rz(-2.2471273) q[3];
sx q[3];
rz(0.47190445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9149949) q[2];
sx q[2];
rz(-1.8151585) q[2];
sx q[2];
rz(2.162852) q[2];
rz(-1.8799051) q[3];
sx q[3];
rz(-2.1224969) q[3];
sx q[3];
rz(-2.252388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9642445) q[0];
sx q[0];
rz(-1.6151936) q[0];
sx q[0];
rz(0.88441315) q[0];
rz(-0.85909596) q[1];
sx q[1];
rz(-2.0180549) q[1];
sx q[1];
rz(-2.9335847) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7743083) q[0];
sx q[0];
rz(-1.2173664) q[0];
sx q[0];
rz(0.53121451) q[0];
rz(-pi) q[1];
x q[1];
rz(1.154083) q[2];
sx q[2];
rz(-2.1325958) q[2];
sx q[2];
rz(2.3903008) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2629548) q[1];
sx q[1];
rz(-2.4391101) q[1];
sx q[1];
rz(1.9874057) q[1];
x q[2];
rz(2.0741803) q[3];
sx q[3];
rz(-1.5690104) q[3];
sx q[3];
rz(-1.747365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.38006833) q[2];
sx q[2];
rz(-0.91975776) q[2];
sx q[2];
rz(0.541614) q[2];
rz(-1.2540865) q[3];
sx q[3];
rz(-1.5897635) q[3];
sx q[3];
rz(2.2327173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8054473) q[0];
sx q[0];
rz(-2.8271524) q[0];
sx q[0];
rz(-2.050515) q[0];
rz(-2.3167141) q[1];
sx q[1];
rz(-2.7179317) q[1];
sx q[1];
rz(1.0728015) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.049041852) q[0];
sx q[0];
rz(-1.9922087) q[0];
sx q[0];
rz(1.7803935) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0457951) q[2];
sx q[2];
rz(-3.1154251) q[2];
sx q[2];
rz(1.2452084) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3531216) q[1];
sx q[1];
rz(-1.4375198) q[1];
sx q[1];
rz(2.8972096) q[1];
x q[2];
rz(2.0951525) q[3];
sx q[3];
rz(-0.45733157) q[3];
sx q[3];
rz(-0.90378896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9020646) q[2];
sx q[2];
rz(-0.54486474) q[2];
sx q[2];
rz(5/(6*pi)) q[2];
rz(-0.15051633) q[3];
sx q[3];
rz(-1.1082606) q[3];
sx q[3];
rz(-0.36293852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6657669) q[0];
sx q[0];
rz(-0.095223991) q[0];
sx q[0];
rz(-1.6640523) q[0];
rz(-2.3382969) q[1];
sx q[1];
rz(-1.3731615) q[1];
sx q[1];
rz(-2.1243748) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8724339) q[0];
sx q[0];
rz(-0.24155051) q[0];
sx q[0];
rz(1.0046602) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3147589) q[2];
sx q[2];
rz(-2.0190329) q[2];
sx q[2];
rz(-2.2729276) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2443631) q[1];
sx q[1];
rz(-1.23037) q[1];
sx q[1];
rz(1.2610408) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3961762) q[3];
sx q[3];
rz(-2.7099797) q[3];
sx q[3];
rz(3.0109757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.21641714) q[2];
sx q[2];
rz(-1.567013) q[2];
sx q[2];
rz(-2.4851921) q[2];
rz(-0.21550719) q[3];
sx q[3];
rz(-1.4114722) q[3];
sx q[3];
rz(1.6925252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.086275252) q[0];
sx q[0];
rz(-1.3441939) q[0];
sx q[0];
rz(-2.2917746) q[0];
rz(2.1814116) q[1];
sx q[1];
rz(-0.4147059) q[1];
sx q[1];
rz(-2.2023831) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.98039) q[0];
sx q[0];
rz(-1.8341258) q[0];
sx q[0];
rz(-0.80391802) q[0];
x q[1];
rz(-1.1562111) q[2];
sx q[2];
rz(-1.5636946) q[2];
sx q[2];
rz(0.19609253) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.43329217) q[1];
sx q[1];
rz(-2.1850249) q[1];
sx q[1];
rz(-1.1105762) q[1];
rz(-pi) q[2];
rz(1.0435095) q[3];
sx q[3];
rz(-1.1498972) q[3];
sx q[3];
rz(2.7038006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.5592929) q[2];
sx q[2];
rz(-1.7377995) q[2];
sx q[2];
rz(0.55029184) q[2];
rz(3.0136717) q[3];
sx q[3];
rz(-3.0030799) q[3];
sx q[3];
rz(0.96854717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48548231) q[0];
sx q[0];
rz(-0.67018296) q[0];
sx q[0];
rz(-0.67076587) q[0];
rz(0.35519629) q[1];
sx q[1];
rz(-1.6658446) q[1];
sx q[1];
rz(2.7459941) q[1];
rz(-0.14329362) q[2];
sx q[2];
rz(-1.431965) q[2];
sx q[2];
rz(1.9964249) q[2];
rz(2.0215423) q[3];
sx q[3];
rz(-1.3629631) q[3];
sx q[3];
rz(-0.25615389) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
