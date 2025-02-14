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
rz(-0.40989447) q[0];
sx q[0];
rz(4.8978187) q[0];
sx q[0];
rz(10.495821) q[0];
rz(1.8154124) q[1];
sx q[1];
rz(2.2161127) q[1];
sx q[1];
rz(9.767017) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15187787) q[0];
sx q[0];
rz(-2.3513455) q[0];
sx q[0];
rz(1.0343767) q[0];
rz(0.81011198) q[2];
sx q[2];
rz(-1.9829921) q[2];
sx q[2];
rz(0.10998943) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6733779) q[1];
sx q[1];
rz(-1.6771375) q[1];
sx q[1];
rz(-2.0061226) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1245021) q[3];
sx q[3];
rz(-1.1747163) q[3];
sx q[3];
rz(-1.210496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5388415) q[2];
sx q[2];
rz(-0.72776908) q[2];
sx q[2];
rz(-1.7068498) q[2];
rz(-2.1753066) q[3];
sx q[3];
rz(-1.1937701) q[3];
sx q[3];
rz(-0.57893354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
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
rz(0.12427881) q[0];
sx q[0];
rz(-0.74324981) q[0];
sx q[0];
rz(3.0857575) q[0];
rz(-0.75885478) q[1];
sx q[1];
rz(-1.1412303) q[1];
sx q[1];
rz(-2.8208044) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0017805) q[0];
sx q[0];
rz(-0.81808144) q[0];
sx q[0];
rz(-1.3670188) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9174395) q[2];
sx q[2];
rz(-0.8220368) q[2];
sx q[2];
rz(1.3393266) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.546678) q[1];
sx q[1];
rz(-1.8045097) q[1];
sx q[1];
rz(-3.1114034) q[1];
x q[2];
rz(1.1029712) q[3];
sx q[3];
rz(-0.89455869) q[3];
sx q[3];
rz(0.53550628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.32030973) q[2];
sx q[2];
rz(-1.8822957) q[2];
sx q[2];
rz(2.8584976) q[2];
rz(-2.8458332) q[3];
sx q[3];
rz(-1.4247954) q[3];
sx q[3];
rz(-1.504771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54842424) q[0];
sx q[0];
rz(-1.230509) q[0];
sx q[0];
rz(-3.044627) q[0];
rz(-1.7123669) q[1];
sx q[1];
rz(-1.2173419) q[1];
sx q[1];
rz(2.7860577) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.36847) q[0];
sx q[0];
rz(-0.63798824) q[0];
sx q[0];
rz(-0.49686749) q[0];
x q[1];
rz(-2.9121551) q[2];
sx q[2];
rz(-2.7167453) q[2];
sx q[2];
rz(-0.72553023) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.83427338) q[1];
sx q[1];
rz(-1.4774228) q[1];
sx q[1];
rz(-1.982467) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.25162752) q[3];
sx q[3];
rz(-1.6109924) q[3];
sx q[3];
rz(-3.0305559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.43704438) q[2];
sx q[2];
rz(-1.9883678) q[2];
sx q[2];
rz(-2.8538749) q[2];
rz(0.25117609) q[3];
sx q[3];
rz(-2.0468678) q[3];
sx q[3];
rz(1.8672966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9525725) q[0];
sx q[0];
rz(-0.59307161) q[0];
sx q[0];
rz(0.75743842) q[0];
rz(-2.3144552) q[1];
sx q[1];
rz(-2.4816781) q[1];
sx q[1];
rz(1.8338902) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8104202) q[0];
sx q[0];
rz(-3.0068827) q[0];
sx q[0];
rz(1.9707546) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8355952) q[2];
sx q[2];
rz(-1.1192516) q[2];
sx q[2];
rz(-1.2689607) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8825619) q[1];
sx q[1];
rz(-2.3271431) q[1];
sx q[1];
rz(2.7689054) q[1];
rz(0.060606011) q[3];
sx q[3];
rz(-1.8148141) q[3];
sx q[3];
rz(0.054494245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7149675) q[2];
sx q[2];
rz(-1.3447309) q[2];
sx q[2];
rz(1.3955383) q[2];
rz(1.3501781) q[3];
sx q[3];
rz(-1.6387286) q[3];
sx q[3];
rz(0.028039886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79513079) q[0];
sx q[0];
rz(-0.53855723) q[0];
sx q[0];
rz(-1.441347) q[0];
rz(-0.92789188) q[1];
sx q[1];
rz(-2.3065232) q[1];
sx q[1];
rz(-0.79383129) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2071811) q[0];
sx q[0];
rz(-2.7156805) q[0];
sx q[0];
rz(1.3176729) q[0];
rz(0.16818856) q[2];
sx q[2];
rz(-1.6861746) q[2];
sx q[2];
rz(-2.4720129) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.87463996) q[1];
sx q[1];
rz(-0.62332223) q[1];
sx q[1];
rz(-1.0250574) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6364274) q[3];
sx q[3];
rz(-1.9117711) q[3];
sx q[3];
rz(-2.7974432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.64357197) q[2];
sx q[2];
rz(-1.6732432) q[2];
sx q[2];
rz(2.7872861) q[2];
rz(-2.0453359) q[3];
sx q[3];
rz(-0.26429629) q[3];
sx q[3];
rz(1.8335584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
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
rz(-0.7588014) q[0];
sx q[0];
rz(-2.379874) q[0];
sx q[0];
rz(2.288901) q[0];
rz(0.26607749) q[1];
sx q[1];
rz(-2.1177025) q[1];
sx q[1];
rz(-1.783225) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5772029) q[0];
sx q[0];
rz(-1.7856585) q[0];
sx q[0];
rz(-0.99265088) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2439315) q[2];
sx q[2];
rz(-1.7836469) q[2];
sx q[2];
rz(1.1291248) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3682421) q[1];
sx q[1];
rz(-2.1424286) q[1];
sx q[1];
rz(-1.9707457) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.27652506) q[3];
sx q[3];
rz(-0.92036906) q[3];
sx q[3];
rz(0.98052935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.69409662) q[2];
sx q[2];
rz(-1.9393549) q[2];
sx q[2];
rz(-0.48074943) q[2];
rz(-0.9350183) q[3];
sx q[3];
rz(-2.1604249) q[3];
sx q[3];
rz(0.41486129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9347436) q[0];
sx q[0];
rz(-2.594279) q[0];
sx q[0];
rz(-0.68688399) q[0];
rz(0.96784776) q[1];
sx q[1];
rz(-1.6136439) q[1];
sx q[1];
rz(2.9639249) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2813276) q[0];
sx q[0];
rz(-1.9138558) q[0];
sx q[0];
rz(2.706091) q[0];
rz(-pi) q[1];
rz(-2.3697621) q[2];
sx q[2];
rz(-0.83862156) q[2];
sx q[2];
rz(3.0385303) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.93157576) q[1];
sx q[1];
rz(-2.615395) q[1];
sx q[1];
rz(-0.99858649) q[1];
x q[2];
rz(-0.92063825) q[3];
sx q[3];
rz(-3.0282109) q[3];
sx q[3];
rz(-0.46550289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.054242) q[2];
sx q[2];
rz(-0.44173104) q[2];
sx q[2];
rz(-2.1596215) q[2];
rz(-2.8864077) q[3];
sx q[3];
rz(-1.1533777) q[3];
sx q[3];
rz(1.0758146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0725237) q[0];
sx q[0];
rz(-0.60966063) q[0];
sx q[0];
rz(-3.0141818) q[0];
rz(-1.6690856) q[1];
sx q[1];
rz(-1.1839048) q[1];
sx q[1];
rz(1.6944108) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5965737) q[0];
sx q[0];
rz(-1.2535639) q[0];
sx q[0];
rz(1.7659643) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9451649) q[2];
sx q[2];
rz(-1.0584694) q[2];
sx q[2];
rz(0.62498876) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.060272) q[1];
sx q[1];
rz(-1.6439207) q[1];
sx q[1];
rz(0.13797073) q[1];
rz(-1.9384222) q[3];
sx q[3];
rz(-0.6163839) q[3];
sx q[3];
rz(-0.65983696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.89333263) q[2];
sx q[2];
rz(-2.1620763) q[2];
sx q[2];
rz(2.7323006) q[2];
rz(-2.4457757) q[3];
sx q[3];
rz(-0.69605637) q[3];
sx q[3];
rz(-2.5223562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61243397) q[0];
sx q[0];
rz(-0.22888628) q[0];
sx q[0];
rz(-2.9869475) q[0];
rz(-1.0908499) q[1];
sx q[1];
rz(-0.88328528) q[1];
sx q[1];
rz(1.326391) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94448254) q[0];
sx q[0];
rz(-1.3061227) q[0];
sx q[0];
rz(2.2933084) q[0];
x q[1];
rz(1.6589734) q[2];
sx q[2];
rz(-0.50082654) q[2];
sx q[2];
rz(-1.0732936) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0834956) q[1];
sx q[1];
rz(-0.82193437) q[1];
sx q[1];
rz(0.067598299) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6821413) q[3];
sx q[3];
rz(-1.3927476) q[3];
sx q[3];
rz(-1.1332172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.078330127) q[2];
sx q[2];
rz(-2.7276954) q[2];
sx q[2];
rz(2.4388893) q[2];
rz(3.0472158) q[3];
sx q[3];
rz(-0.97797147) q[3];
sx q[3];
rz(0.92488658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6150045) q[0];
sx q[0];
rz(-1.0300535) q[0];
sx q[0];
rz(-2.7870542) q[0];
rz(-1.7189369) q[1];
sx q[1];
rz(-1.4965897) q[1];
sx q[1];
rz(2.6197701) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3799135) q[0];
sx q[0];
rz(-1.2186134) q[0];
sx q[0];
rz(-0.16009472) q[0];
x q[1];
rz(-1.8759236) q[2];
sx q[2];
rz(-1.4386476) q[2];
sx q[2];
rz(-0.73763359) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3067535) q[1];
sx q[1];
rz(-2.2242079) q[1];
sx q[1];
rz(0.70717137) q[1];
rz(-pi) q[2];
x q[2];
rz(0.75504889) q[3];
sx q[3];
rz(-2.2590989) q[3];
sx q[3];
rz(1.2699708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3598513) q[2];
sx q[2];
rz(-1.1726817) q[2];
sx q[2];
rz(0.68217984) q[2];
rz(-1.773268) q[3];
sx q[3];
rz(-1.6006288) q[3];
sx q[3];
rz(-2.1413596) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2405887) q[0];
sx q[0];
rz(-0.95745845) q[0];
sx q[0];
rz(-0.29314713) q[0];
rz(0.82943574) q[1];
sx q[1];
rz(-1.4301626) q[1];
sx q[1];
rz(1.593874) q[1];
rz(2.9113795) q[2];
sx q[2];
rz(-1.3166192) q[2];
sx q[2];
rz(-0.26659378) q[2];
rz(1.6080018) q[3];
sx q[3];
rz(-1.6757575) q[3];
sx q[3];
rz(2.5057275) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
