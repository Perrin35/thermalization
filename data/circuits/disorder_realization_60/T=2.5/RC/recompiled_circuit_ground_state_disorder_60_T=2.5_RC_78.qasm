OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-3.056501) q[0];
sx q[0];
rz(-2.8488475) q[0];
sx q[0];
rz(-2.2226287) q[0];
rz(2.7594944) q[1];
sx q[1];
rz(-0.16799071) q[1];
sx q[1];
rz(-1.1408495) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1534682) q[0];
sx q[0];
rz(-1.690295) q[0];
sx q[0];
rz(-1.7480127) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4130526) q[2];
sx q[2];
rz(-2.72914) q[2];
sx q[2];
rz(2.0714687) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7431822) q[1];
sx q[1];
rz(-0.40603058) q[1];
sx q[1];
rz(2.7999785) q[1];
x q[2];
rz(-1.857418) q[3];
sx q[3];
rz(-0.45551963) q[3];
sx q[3];
rz(2.974433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4484278) q[2];
sx q[2];
rz(-2.0824671) q[2];
sx q[2];
rz(1.9654174) q[2];
rz(0.042393427) q[3];
sx q[3];
rz(-1.8040801) q[3];
sx q[3];
rz(-0.5347518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(2.6723044) q[0];
sx q[0];
rz(-0.35457087) q[0];
sx q[0];
rz(0.18445036) q[0];
rz(3.0448659) q[1];
sx q[1];
rz(-1.5238785) q[1];
sx q[1];
rz(-1.2299889) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1095268) q[0];
sx q[0];
rz(-1.7735039) q[0];
sx q[0];
rz(-0.37287449) q[0];
rz(2.582091) q[2];
sx q[2];
rz(-1.6116665) q[2];
sx q[2];
rz(-2.3377583) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.36791554) q[1];
sx q[1];
rz(-1.8080447) q[1];
sx q[1];
rz(1.3648454) q[1];
rz(-pi) q[2];
rz(0.038944728) q[3];
sx q[3];
rz(-1.9488261) q[3];
sx q[3];
rz(1.8824487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7868598) q[2];
sx q[2];
rz(-2.730098) q[2];
sx q[2];
rz(0.011431781) q[2];
rz(-1.8276021) q[3];
sx q[3];
rz(-1.7766137) q[3];
sx q[3];
rz(0.7029117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3979724) q[0];
sx q[0];
rz(-0.13298661) q[0];
sx q[0];
rz(0.61070329) q[0];
rz(-0.01344219) q[1];
sx q[1];
rz(-1.2862658) q[1];
sx q[1];
rz(0.59649831) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0749576) q[0];
sx q[0];
rz(-0.63371113) q[0];
sx q[0];
rz(1.2873961) q[0];
rz(0.99278583) q[2];
sx q[2];
rz(-1.3137378) q[2];
sx q[2];
rz(-1.7456733) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.33112673) q[1];
sx q[1];
rz(-2.2743723) q[1];
sx q[1];
rz(-1.1570279) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5174808) q[3];
sx q[3];
rz(-0.55453306) q[3];
sx q[3];
rz(2.7729101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.67768031) q[2];
sx q[2];
rz(-1.3119421) q[2];
sx q[2];
rz(0.00027351969) q[2];
rz(1.6655946) q[3];
sx q[3];
rz(-0.6690343) q[3];
sx q[3];
rz(0.10703787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45604712) q[0];
sx q[0];
rz(-0.16655971) q[0];
sx q[0];
rz(-1.1628994) q[0];
rz(-2.1641425) q[1];
sx q[1];
rz(-1.5403427) q[1];
sx q[1];
rz(-1.5553442) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9396203) q[0];
sx q[0];
rz(-0.012379025) q[0];
sx q[0];
rz(-0.89039652) q[0];
rz(-pi) q[1];
rz(1.5260209) q[2];
sx q[2];
rz(-0.78592671) q[2];
sx q[2];
rz(1.3040552) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.11862893) q[1];
sx q[1];
rz(-1.5852155) q[1];
sx q[1];
rz(-0.49683907) q[1];
rz(-2.0647207) q[3];
sx q[3];
rz(-1.3528429) q[3];
sx q[3];
rz(-2.2214784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1996475) q[2];
sx q[2];
rz(-0.5117828) q[2];
sx q[2];
rz(-0.23435782) q[2];
rz(0.50218454) q[3];
sx q[3];
rz(-1.0415223) q[3];
sx q[3];
rz(-0.65619367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31768826) q[0];
sx q[0];
rz(-0.81517878) q[0];
sx q[0];
rz(1.748388) q[0];
rz(2.4488917) q[1];
sx q[1];
rz(-1.0857948) q[1];
sx q[1];
rz(2.0511973) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4596953) q[0];
sx q[0];
rz(-2.5715264) q[0];
sx q[0];
rz(1.9697813) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2545891) q[2];
sx q[2];
rz(-0.74216026) q[2];
sx q[2];
rz(0.61486828) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7140035) q[1];
sx q[1];
rz(-1.98381) q[1];
sx q[1];
rz(-0.16431134) q[1];
rz(-pi) q[2];
rz(1.0308068) q[3];
sx q[3];
rz(-1.705369) q[3];
sx q[3];
rz(-2.7887044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.580487) q[2];
sx q[2];
rz(-1.9876391) q[2];
sx q[2];
rz(-0.53606501) q[2];
rz(-2.8481893) q[3];
sx q[3];
rz(-1.2111726) q[3];
sx q[3];
rz(2.6611879) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2368161) q[0];
sx q[0];
rz(-0.95229709) q[0];
sx q[0];
rz(3.0425518) q[0];
rz(-2.1095443) q[1];
sx q[1];
rz(-2.2910304) q[1];
sx q[1];
rz(1.438407) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0625298) q[0];
sx q[0];
rz(-2.2465116) q[0];
sx q[0];
rz(2.8131301) q[0];
rz(-2.8581335) q[2];
sx q[2];
rz(-1.6416896) q[2];
sx q[2];
rz(2.9170582) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9248396) q[1];
sx q[1];
rz(-0.79063639) q[1];
sx q[1];
rz(3.1086712) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.766633) q[3];
sx q[3];
rz(-0.31878933) q[3];
sx q[3];
rz(0.44170435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9859887) q[2];
sx q[2];
rz(-2.6376548) q[2];
sx q[2];
rz(-0.8209374) q[2];
rz(-1.692903) q[3];
sx q[3];
rz(-2.4235348) q[3];
sx q[3];
rz(-1.6528355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2445225) q[0];
sx q[0];
rz(-0.91461602) q[0];
sx q[0];
rz(-2.6389417) q[0];
rz(1.2333168) q[1];
sx q[1];
rz(-2.2063467) q[1];
sx q[1];
rz(-2.1615692) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5222447) q[0];
sx q[0];
rz(-2.7947593) q[0];
sx q[0];
rz(0.808544) q[0];
x q[1];
rz(-1.8623615) q[2];
sx q[2];
rz(-0.81613651) q[2];
sx q[2];
rz(-1.3929491) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7240643) q[1];
sx q[1];
rz(-1.8019946) q[1];
sx q[1];
rz(2.3910012) q[1];
rz(-1.8007038) q[3];
sx q[3];
rz(-0.74384825) q[3];
sx q[3];
rz(-0.4901674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9194453) q[2];
sx q[2];
rz(-1.8541226) q[2];
sx q[2];
rz(-1.7670828) q[2];
rz(-1.3162656) q[3];
sx q[3];
rz(-2.2856183) q[3];
sx q[3];
rz(2.159806) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.114349) q[0];
sx q[0];
rz(-2.7482432) q[0];
sx q[0];
rz(2.223176) q[0];
rz(-2.9367327) q[1];
sx q[1];
rz(-1.058895) q[1];
sx q[1];
rz(-1.4835666) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99824725) q[0];
sx q[0];
rz(-1.424331) q[0];
sx q[0];
rz(-1.2781758) q[0];
rz(-pi) q[1];
rz(1.7006364) q[2];
sx q[2];
rz(-2.3774638) q[2];
sx q[2];
rz(1.6418599) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.10459722) q[1];
sx q[1];
rz(-0.87338398) q[1];
sx q[1];
rz(-0.11499494) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0851787) q[3];
sx q[3];
rz(-2.4216363) q[3];
sx q[3];
rz(2.3145793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1398937) q[2];
sx q[2];
rz(-0.49109444) q[2];
sx q[2];
rz(1.3341382) q[2];
rz(-2.259518) q[3];
sx q[3];
rz(-2.4460654) q[3];
sx q[3];
rz(1.2923366) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.081721574) q[0];
sx q[0];
rz(-2.7099755) q[0];
sx q[0];
rz(0.01817848) q[0];
rz(0.40953088) q[1];
sx q[1];
rz(-1.7117701) q[1];
sx q[1];
rz(3.0407564) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6438599) q[0];
sx q[0];
rz(-2.5063305) q[0];
sx q[0];
rz(0.058202581) q[0];
x q[1];
rz(-1.6214431) q[2];
sx q[2];
rz(-1.9910004) q[2];
sx q[2];
rz(0.24383185) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3506895) q[1];
sx q[1];
rz(-1.8115582) q[1];
sx q[1];
rz(1.945182) q[1];
rz(-pi) q[2];
rz(3.0206817) q[3];
sx q[3];
rz(-1.8163346) q[3];
sx q[3];
rz(0.80805486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.33710256) q[2];
sx q[2];
rz(-3.0323995) q[2];
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
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(0.2357904) q[0];
sx q[0];
rz(-1.6658655) q[0];
sx q[0];
rz(-2.8210848) q[0];
rz(-0.39101741) q[1];
sx q[1];
rz(-1.1415488) q[1];
sx q[1];
rz(1.5919707) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79637291) q[0];
sx q[0];
rz(-1.380305) q[0];
sx q[0];
rz(-0.20261554) q[0];
x q[1];
rz(-2.6274849) q[2];
sx q[2];
rz(-1.7127345) q[2];
sx q[2];
rz(-0.98279233) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1913073) q[1];
sx q[1];
rz(-0.37883329) q[1];
sx q[1];
rz(-1.8089507) q[1];
x q[2];
rz(2.7192468) q[3];
sx q[3];
rz(-1.5639439) q[3];
sx q[3];
rz(0.50240483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.090342) q[2];
sx q[2];
rz(-2.0917459) q[2];
sx q[2];
rz(2.7412097) q[2];
rz(2.9054437) q[3];
sx q[3];
rz(-3.0381687) q[3];
sx q[3];
rz(2.1307438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26265963) q[0];
sx q[0];
rz(-2.6317609) q[0];
sx q[0];
rz(1.2608933) q[0];
rz(2.6564468) q[1];
sx q[1];
rz(-0.79882516) q[1];
sx q[1];
rz(-0.47732236) q[1];
rz(-3.015949) q[2];
sx q[2];
rz(-1.3640654) q[2];
sx q[2];
rz(-2.9403093) q[2];
rz(1.7988206) q[3];
sx q[3];
rz(-1.5088456) q[3];
sx q[3];
rz(0.07901293) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
