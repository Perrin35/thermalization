OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.52842927) q[0];
sx q[0];
rz(-1.0597205) q[0];
sx q[0];
rz(-2.4106195) q[0];
rz(-1.5001186) q[1];
sx q[1];
rz(-2.1067696) q[1];
sx q[1];
rz(-2.1980481) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65981728) q[0];
sx q[0];
rz(-3.0660015) q[0];
sx q[0];
rz(-2.6600921) q[0];
x q[1];
rz(1.4859096) q[2];
sx q[2];
rz(-0.92157084) q[2];
sx q[2];
rz(-0.50253403) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.01602068) q[1];
sx q[1];
rz(-1.9106094) q[1];
sx q[1];
rz(1.1556975) q[1];
x q[2];
rz(-2.5991873) q[3];
sx q[3];
rz(-2.7895658) q[3];
sx q[3];
rz(2.5792518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8865108) q[2];
sx q[2];
rz(-1.7604897) q[2];
sx q[2];
rz(1.8908267) q[2];
rz(1.7154153) q[3];
sx q[3];
rz(-0.91606796) q[3];
sx q[3];
rz(-0.9799408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.132906) q[0];
sx q[0];
rz(-2.0693021) q[0];
sx q[0];
rz(-2.5426478) q[0];
rz(-1.8006181) q[1];
sx q[1];
rz(-0.95021617) q[1];
sx q[1];
rz(-0.96639955) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37106284) q[0];
sx q[0];
rz(-2.3822228) q[0];
sx q[0];
rz(0.66803996) q[0];
rz(-1.5210549) q[2];
sx q[2];
rz(-0.94919862) q[2];
sx q[2];
rz(3.0034686) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4557138) q[1];
sx q[1];
rz(-3.0027632) q[1];
sx q[1];
rz(-0.60771897) q[1];
rz(-pi) q[2];
x q[2];
rz(0.74921272) q[3];
sx q[3];
rz(-1.2470761) q[3];
sx q[3];
rz(2.2183228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0559343) q[2];
sx q[2];
rz(-0.81515437) q[2];
sx q[2];
rz(-0.30109626) q[2];
rz(-1.9484693) q[3];
sx q[3];
rz(-1.5501225) q[3];
sx q[3];
rz(-1.1863856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0369204) q[0];
sx q[0];
rz(-1.6372697) q[0];
sx q[0];
rz(-1.954129) q[0];
rz(-1.9056412) q[1];
sx q[1];
rz(-1.0373479) q[1];
sx q[1];
rz(-1.8240066) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7147489) q[0];
sx q[0];
rz(-2.2509529) q[0];
sx q[0];
rz(-1.2965409) q[0];
x q[1];
rz(-0.74430978) q[2];
sx q[2];
rz(-2.8544606) q[2];
sx q[2];
rz(-0.14936514) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6685851) q[1];
sx q[1];
rz(-2.7088532) q[1];
sx q[1];
rz(2.8981478) q[1];
x q[2];
rz(0.1173238) q[3];
sx q[3];
rz(-2.085272) q[3];
sx q[3];
rz(2.8770212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0111771) q[2];
sx q[2];
rz(-1.3935564) q[2];
sx q[2];
rz(-0.2066361) q[2];
rz(-2.4335499) q[3];
sx q[3];
rz(-2.9338624) q[3];
sx q[3];
rz(-0.89282435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77465039) q[0];
sx q[0];
rz(-0.21454021) q[0];
sx q[0];
rz(2.2553717) q[0];
rz(-1.0097424) q[1];
sx q[1];
rz(-0.90615288) q[1];
sx q[1];
rz(1.9151691) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7851631) q[0];
sx q[0];
rz(-2.4117081) q[0];
sx q[0];
rz(-1.6701783) q[0];
x q[1];
rz(-0.88422758) q[2];
sx q[2];
rz(-1.9366169) q[2];
sx q[2];
rz(-0.58931749) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9981421) q[1];
sx q[1];
rz(-1.3610024) q[1];
sx q[1];
rz(1.7721121) q[1];
x q[2];
rz(1.1981443) q[3];
sx q[3];
rz(-2.028392) q[3];
sx q[3];
rz(0.83329337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4975138) q[2];
sx q[2];
rz(-1.8277233) q[2];
sx q[2];
rz(-0.99299661) q[2];
rz(-1.8289061) q[3];
sx q[3];
rz(-1.0220746) q[3];
sx q[3];
rz(-0.48373568) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.118367) q[0];
sx q[0];
rz(-0.3148196) q[0];
sx q[0];
rz(-1.9556048) q[0];
rz(1.7182619) q[1];
sx q[1];
rz(-1.490373) q[1];
sx q[1];
rz(-0.58247724) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0156508) q[0];
sx q[0];
rz(-1.6424718) q[0];
sx q[0];
rz(-2.1732251) q[0];
rz(-pi) q[1];
rz(0.8930348) q[2];
sx q[2];
rz(-1.4544011) q[2];
sx q[2];
rz(1.900577) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7354048) q[1];
sx q[1];
rz(-1.5516073) q[1];
sx q[1];
rz(-0.71449844) q[1];
rz(1.6378239) q[3];
sx q[3];
rz(-1.9762632) q[3];
sx q[3];
rz(-2.1342579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4867268) q[2];
sx q[2];
rz(-1.5348237) q[2];
sx q[2];
rz(3.0598818) q[2];
rz(2.667526) q[3];
sx q[3];
rz(-1.2636377) q[3];
sx q[3];
rz(-1.6430395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24494568) q[0];
sx q[0];
rz(-2.0136254) q[0];
sx q[0];
rz(3.1337877) q[0];
rz(-1.4004978) q[1];
sx q[1];
rz(-0.85406071) q[1];
sx q[1];
rz(-2.0369464) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68708006) q[0];
sx q[0];
rz(-0.91327635) q[0];
sx q[0];
rz(-1.7772654) q[0];
rz(-pi) q[1];
rz(-1.3278264) q[2];
sx q[2];
rz(-0.94488482) q[2];
sx q[2];
rz(-1.5205795) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1659516) q[1];
sx q[1];
rz(-0.8756606) q[1];
sx q[1];
rz(2.402311) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1052746) q[3];
sx q[3];
rz(-1.7200617) q[3];
sx q[3];
rz(-0.6061337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2509987) q[2];
sx q[2];
rz(-0.72967523) q[2];
sx q[2];
rz(2.5863623) q[2];
rz(-0.17272078) q[3];
sx q[3];
rz(-1.828086) q[3];
sx q[3];
rz(-1.593332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5884488) q[0];
sx q[0];
rz(-0.48184904) q[0];
sx q[0];
rz(3.0798262) q[0];
rz(-0.24208367) q[1];
sx q[1];
rz(-0.37529072) q[1];
sx q[1];
rz(-1.1118836) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0985078) q[0];
sx q[0];
rz(-2.3102009) q[0];
sx q[0];
rz(2.1466473) q[0];
rz(3.1277666) q[2];
sx q[2];
rz(-1.4931884) q[2];
sx q[2];
rz(1.3071878) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.71373938) q[1];
sx q[1];
rz(-2.3728328) q[1];
sx q[1];
rz(-2.7379235) q[1];
rz(-pi) q[2];
rz(-0.32825177) q[3];
sx q[3];
rz(-1.8040856) q[3];
sx q[3];
rz(2.1094028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.8093439) q[2];
sx q[2];
rz(-0.76433864) q[2];
sx q[2];
rz(-0.81364441) q[2];
rz(-1.404473) q[3];
sx q[3];
rz(-2.9128894) q[3];
sx q[3];
rz(-2.5261734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62548816) q[0];
sx q[0];
rz(-1.7913211) q[0];
sx q[0];
rz(-1.7161436) q[0];
rz(1.5215993) q[1];
sx q[1];
rz(-0.75192538) q[1];
sx q[1];
rz(-0.63751784) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91598788) q[0];
sx q[0];
rz(-0.84945744) q[0];
sx q[0];
rz(2.3774873) q[0];
rz(-0.72458467) q[2];
sx q[2];
rz(-1.2983592) q[2];
sx q[2];
rz(0.1833293) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.058640826) q[1];
sx q[1];
rz(-2.4070027) q[1];
sx q[1];
rz(-2.525108) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7082801) q[3];
sx q[3];
rz(-2.1372037) q[3];
sx q[3];
rz(-2.8979104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1311538) q[2];
sx q[2];
rz(-0.70019478) q[2];
sx q[2];
rz(-2.0054224) q[2];
rz(1.6561967) q[3];
sx q[3];
rz(-2.6172726) q[3];
sx q[3];
rz(-1.9320528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5159601) q[0];
sx q[0];
rz(-0.78173286) q[0];
sx q[0];
rz(-0.67614722) q[0];
rz(0.82538429) q[1];
sx q[1];
rz(-2.717658) q[1];
sx q[1];
rz(2.1264145) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0494173) q[0];
sx q[0];
rz(-1.1560688) q[0];
sx q[0];
rz(-1.7134922) q[0];
x q[1];
rz(-1.4170309) q[2];
sx q[2];
rz(-0.92369881) q[2];
sx q[2];
rz(0.53714067) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.79267348) q[1];
sx q[1];
rz(-1.9209849) q[1];
sx q[1];
rz(0.37383553) q[1];
rz(-pi) q[2];
rz(-2.7916662) q[3];
sx q[3];
rz(-2.3438128) q[3];
sx q[3];
rz(1.6328904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4201346) q[2];
sx q[2];
rz(-2.9253503) q[2];
sx q[2];
rz(-0.96735111) q[2];
rz(-1.5445276) q[3];
sx q[3];
rz(-1.8287851) q[3];
sx q[3];
rz(-2.7887204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5230781) q[0];
sx q[0];
rz(-1.5788364) q[0];
sx q[0];
rz(3.0850947) q[0];
rz(2.0691195) q[1];
sx q[1];
rz(-1.871855) q[1];
sx q[1];
rz(-1.7369695) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0427986) q[0];
sx q[0];
rz(-0.3846752) q[0];
sx q[0];
rz(1.2205475) q[0];
x q[1];
rz(1.6925473) q[2];
sx q[2];
rz(-1.2687614) q[2];
sx q[2];
rz(-1.2388602) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.9359293) q[1];
sx q[1];
rz(-0.91055369) q[1];
sx q[1];
rz(-0.6026938) q[1];
rz(-pi) q[2];
rz(0.057617188) q[3];
sx q[3];
rz(-0.15078292) q[3];
sx q[3];
rz(-2.5654716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.29356062) q[2];
sx q[2];
rz(-0.75389391) q[2];
sx q[2];
rz(0.30612293) q[2];
rz(2.7434769) q[3];
sx q[3];
rz(-1.6653776) q[3];
sx q[3];
rz(-1.0242296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8067779) q[0];
sx q[0];
rz(-1.0659185) q[0];
sx q[0];
rz(0.48068) q[0];
rz(1.1595935) q[1];
sx q[1];
rz(-2.1203142) q[1];
sx q[1];
rz(0.28657985) q[1];
rz(2.878306) q[2];
sx q[2];
rz(-1.8594212) q[2];
sx q[2];
rz(-2.6066305) q[2];
rz(1.0767827) q[3];
sx q[3];
rz(-1.8643338) q[3];
sx q[3];
rz(-2.1716933) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];