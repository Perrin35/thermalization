OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-3.1104133) q[0];
sx q[0];
rz(-2.2194982) q[0];
sx q[0];
rz(-2.3715012) q[0];
rz(-2.4251179) q[1];
sx q[1];
rz(-0.94243503) q[1];
sx q[1];
rz(2.4997349) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9489884) q[0];
sx q[0];
rz(-1.3697213) q[0];
sx q[0];
rz(-1.4549535) q[0];
rz(-0.0014512295) q[2];
sx q[2];
rz(-1.5696758) q[2];
sx q[2];
rz(1.4933153) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.042965021) q[1];
sx q[1];
rz(-0.39448276) q[1];
sx q[1];
rz(-0.41542713) q[1];
rz(-pi) q[2];
rz(-2.3083616) q[3];
sx q[3];
rz(-1.1166443) q[3];
sx q[3];
rz(1.3564701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1610819) q[2];
sx q[2];
rz(-2.0445721) q[2];
sx q[2];
rz(2.3316627) q[2];
rz(-0.03446456) q[3];
sx q[3];
rz(-0.66027111) q[3];
sx q[3];
rz(0.1035498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.506839) q[0];
sx q[0];
rz(-0.50951183) q[0];
sx q[0];
rz(-0.65938812) q[0];
rz(1.7022645) q[1];
sx q[1];
rz(-1.5123475) q[1];
sx q[1];
rz(0.48092458) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60431193) q[0];
sx q[0];
rz(-1.3635646) q[0];
sx q[0];
rz(-1.1234493) q[0];
x q[1];
rz(0.76586313) q[2];
sx q[2];
rz(-0.89855902) q[2];
sx q[2];
rz(1.5352676) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9240504) q[1];
sx q[1];
rz(-0.21044193) q[1];
sx q[1];
rz(-1.4242875) q[1];
x q[2];
rz(-2.4830803) q[3];
sx q[3];
rz(-2.4511946) q[3];
sx q[3];
rz(1.8883324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3769569) q[2];
sx q[2];
rz(-0.83676338) q[2];
sx q[2];
rz(2.5118206) q[2];
rz(1.9624286) q[3];
sx q[3];
rz(-0.7032913) q[3];
sx q[3];
rz(-1.111697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69029194) q[0];
sx q[0];
rz(-3.1307104) q[0];
sx q[0];
rz(0.96845281) q[0];
rz(-0.14006607) q[1];
sx q[1];
rz(-1.7882971) q[1];
sx q[1];
rz(0.5828988) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4752858) q[0];
sx q[0];
rz(-2.9687442) q[0];
sx q[0];
rz(-2.8768507) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.23704657) q[2];
sx q[2];
rz(-2.2893527) q[2];
sx q[2];
rz(-1.3952554) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0839094) q[1];
sx q[1];
rz(-2.1321802) q[1];
sx q[1];
rz(2.3777005) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2517334) q[3];
sx q[3];
rz(-1.8382065) q[3];
sx q[3];
rz(1.2280994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.53050238) q[2];
sx q[2];
rz(-0.077433057) q[2];
sx q[2];
rz(-2.9275242) q[2];
rz(-0.30238447) q[3];
sx q[3];
rz(-0.74367911) q[3];
sx q[3];
rz(-2.7157057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18836235) q[0];
sx q[0];
rz(-2.132405) q[0];
sx q[0];
rz(1.0444214) q[0];
rz(-0.87772477) q[1];
sx q[1];
rz(-1.5526155) q[1];
sx q[1];
rz(-2.9728319) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.343086) q[0];
sx q[0];
rz(-1.0765392) q[0];
sx q[0];
rz(0.25747184) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8977106) q[2];
sx q[2];
rz(-1.3762646) q[2];
sx q[2];
rz(3.1049797) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5280594) q[1];
sx q[1];
rz(-0.643104) q[1];
sx q[1];
rz(2.7763441) q[1];
x q[2];
rz(3.0996347) q[3];
sx q[3];
rz(-0.86554147) q[3];
sx q[3];
rz(0.77159568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3622482) q[2];
sx q[2];
rz(-1.8579973) q[2];
sx q[2];
rz(1.0603504) q[2];
rz(2.849071) q[3];
sx q[3];
rz(-0.7258324) q[3];
sx q[3];
rz(0.084107548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58920687) q[0];
sx q[0];
rz(-0.64738208) q[0];
sx q[0];
rz(0.86828434) q[0];
rz(-2.5153416) q[1];
sx q[1];
rz(-1.3887082) q[1];
sx q[1];
rz(-1.3583604) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.043644301) q[0];
sx q[0];
rz(-1.4820423) q[0];
sx q[0];
rz(0.28907816) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7600438) q[2];
sx q[2];
rz(-1.7722436) q[2];
sx q[2];
rz(0.82009456) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7029801) q[1];
sx q[1];
rz(-1.7963406) q[1];
sx q[1];
rz(1.3425407) q[1];
rz(-pi) q[2];
x q[2];
rz(0.56760889) q[3];
sx q[3];
rz(-0.72101147) q[3];
sx q[3];
rz(-0.22338671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9354349) q[2];
sx q[2];
rz(-2.3249966) q[2];
sx q[2];
rz(2.6182776) q[2];
rz(2.7715136) q[3];
sx q[3];
rz(-2.3612634) q[3];
sx q[3];
rz(-1.3453329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10941457) q[0];
sx q[0];
rz(-2.9258756) q[0];
sx q[0];
rz(0.35414645) q[0];
rz(2.1996563) q[1];
sx q[1];
rz(-1.6862005) q[1];
sx q[1];
rz(1.8249493) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9958651) q[0];
sx q[0];
rz(-2.6168129) q[0];
sx q[0];
rz(0.20916931) q[0];
rz(-2.5923827) q[2];
sx q[2];
rz(-1.5592279) q[2];
sx q[2];
rz(-0.14278296) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1209379) q[1];
sx q[1];
rz(-1.6319425) q[1];
sx q[1];
rz(-0.17423363) q[1];
rz(-pi) q[2];
rz(0.64423465) q[3];
sx q[3];
rz(-1.7984063) q[3];
sx q[3];
rz(1.5952974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8159304) q[2];
sx q[2];
rz(-2.7553835) q[2];
sx q[2];
rz(-0.33829921) q[2];
rz(2.6650688) q[3];
sx q[3];
rz(-0.75348133) q[3];
sx q[3];
rz(0.34059718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7020096) q[0];
sx q[0];
rz(-1.5583353) q[0];
sx q[0];
rz(-2.6038792) q[0];
rz(-1.4078377) q[1];
sx q[1];
rz(-2.6815963) q[1];
sx q[1];
rz(2.5163311) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38442507) q[0];
sx q[0];
rz(-2.129607) q[0];
sx q[0];
rz(-1.1435056) q[0];
rz(-pi) q[1];
rz(-2.602732) q[2];
sx q[2];
rz(-1.5075353) q[2];
sx q[2];
rz(-0.45743313) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(9*pi/10) q[1];
sx q[1];
rz(-2.4971967) q[1];
sx q[1];
rz(-0.36233904) q[1];
rz(1.2042562) q[3];
sx q[3];
rz(-0.82895422) q[3];
sx q[3];
rz(-2.3693565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9739146) q[2];
sx q[2];
rz(-2.6748071) q[2];
sx q[2];
rz(-1.3803049) q[2];
rz(-2.7050833) q[3];
sx q[3];
rz(-2.1197539) q[3];
sx q[3];
rz(-2.4280587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
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
rz(-1.9647144) q[0];
sx q[0];
rz(-2.5202993) q[0];
sx q[0];
rz(-0.51625133) q[0];
rz(-0.77975726) q[1];
sx q[1];
rz(-2.1596491) q[1];
sx q[1];
rz(2.0947184) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31016997) q[0];
sx q[0];
rz(-0.3331694) q[0];
sx q[0];
rz(2.6875671) q[0];
x q[1];
rz(-0.70430906) q[2];
sx q[2];
rz(-0.71137911) q[2];
sx q[2];
rz(-0.34397438) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.88491625) q[1];
sx q[1];
rz(-0.11029989) q[1];
sx q[1];
rz(2.1259456) q[1];
rz(-pi) q[2];
rz(1.9674928) q[3];
sx q[3];
rz(-1.124212) q[3];
sx q[3];
rz(0.74029884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9376935) q[2];
sx q[2];
rz(-2.9475309) q[2];
sx q[2];
rz(0.95721179) q[2];
rz(-0.29414487) q[3];
sx q[3];
rz(-1.1295986) q[3];
sx q[3];
rz(-0.2778151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99943632) q[0];
sx q[0];
rz(-0.23040982) q[0];
sx q[0];
rz(-2.9545422) q[0];
rz(-2.710178) q[1];
sx q[1];
rz(-2.7038733) q[1];
sx q[1];
rz(1.4923219) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.084928345) q[0];
sx q[0];
rz(-1.5001552) q[0];
sx q[0];
rz(-1.4991374) q[0];
x q[1];
rz(-2.951118) q[2];
sx q[2];
rz(-1.2419133) q[2];
sx q[2];
rz(2.637907) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.75892657) q[1];
sx q[1];
rz(-1.1622283) q[1];
sx q[1];
rz(0.17951413) q[1];
rz(-pi) q[2];
rz(-3.047916) q[3];
sx q[3];
rz(-1.4149349) q[3];
sx q[3];
rz(0.81934281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.46818587) q[2];
sx q[2];
rz(-2.2298721) q[2];
sx q[2];
rz(0.98463303) q[2];
rz(-2.6268688) q[3];
sx q[3];
rz(-0.52557164) q[3];
sx q[3];
rz(1.2334067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24432261) q[0];
sx q[0];
rz(-1.6615302) q[0];
sx q[0];
rz(-2.3829714) q[0];
rz(1.1933391) q[1];
sx q[1];
rz(-1.1921644) q[1];
sx q[1];
rz(-1.4512482) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.25989) q[0];
sx q[0];
rz(-2.8326747) q[0];
sx q[0];
rz(2.3007459) q[0];
rz(-pi) q[1];
rz(-1.3682215) q[2];
sx q[2];
rz(-1.5294917) q[2];
sx q[2];
rz(-1.5243204) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7566061) q[1];
sx q[1];
rz(-0.95545095) q[1];
sx q[1];
rz(0.032757515) q[1];
rz(0.32231195) q[3];
sx q[3];
rz(-1.9602284) q[3];
sx q[3];
rz(-3.1403613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.54185581) q[2];
sx q[2];
rz(-0.92076045) q[2];
sx q[2];
rz(0.80121458) q[2];
rz(-0.93252212) q[3];
sx q[3];
rz(-1.9263809) q[3];
sx q[3];
rz(2.7264989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5781317) q[0];
sx q[0];
rz(-1.7628071) q[0];
sx q[0];
rz(2.5264869) q[0];
rz(-2.8201132) q[1];
sx q[1];
rz(-2.165806) q[1];
sx q[1];
rz(1.4478366) q[1];
rz(1.6952487) q[2];
sx q[2];
rz(-1.7348518) q[2];
sx q[2];
rz(3.0748925) q[2];
rz(1.8516171) q[3];
sx q[3];
rz(-1.7924037) q[3];
sx q[3];
rz(-2.0879346) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
