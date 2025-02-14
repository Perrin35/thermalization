OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.428838) q[0];
sx q[0];
rz(-1.5313671) q[0];
sx q[0];
rz(2.03696) q[0];
rz(-1.8072577) q[1];
sx q[1];
rz(-0.72307888) q[1];
sx q[1];
rz(-1.877797) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1570779) q[0];
sx q[0];
rz(-0.70467585) q[0];
sx q[0];
rz(1.9542171) q[0];
rz(-pi) q[1];
x q[1];
rz(0.96226389) q[2];
sx q[2];
rz(-0.66676408) q[2];
sx q[2];
rz(0.61362574) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.25322029) q[1];
sx q[1];
rz(-1.651763) q[1];
sx q[1];
rz(1.8964279) q[1];
x q[2];
rz(-2.0588082) q[3];
sx q[3];
rz(-0.48824874) q[3];
sx q[3];
rz(-1.6475208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9870712) q[2];
sx q[2];
rz(-1.0769341) q[2];
sx q[2];
rz(1.3686352) q[2];
rz(-0.23855071) q[3];
sx q[3];
rz(-2.7261901) q[3];
sx q[3];
rz(0.11428782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40047613) q[0];
sx q[0];
rz(-2.0023161) q[0];
sx q[0];
rz(1.1503295) q[0];
rz(3.053983) q[1];
sx q[1];
rz(-1.3299512) q[1];
sx q[1];
rz(1.5706583) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8920736) q[0];
sx q[0];
rz(-0.49580916) q[0];
sx q[0];
rz(2.9167988) q[0];
rz(-pi) q[1];
x q[1];
rz(0.3377487) q[2];
sx q[2];
rz(-2.9523473) q[2];
sx q[2];
rz(2.7836329) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1024061) q[1];
sx q[1];
rz(-2.3433609) q[1];
sx q[1];
rz(0.95603966) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6704795) q[3];
sx q[3];
rz(-2.2918252) q[3];
sx q[3];
rz(1.6278933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.34438434) q[2];
sx q[2];
rz(-2.4281561) q[2];
sx q[2];
rz(0.1864645) q[2];
rz(0.79948419) q[3];
sx q[3];
rz(-1.5954285) q[3];
sx q[3];
rz(2.1605087) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8568273) q[0];
sx q[0];
rz(-1.3815877) q[0];
sx q[0];
rz(-3.0322266) q[0];
rz(-2.5378387) q[1];
sx q[1];
rz(-1.8021288) q[1];
sx q[1];
rz(1.0341136) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84635705) q[0];
sx q[0];
rz(-0.24938008) q[0];
sx q[0];
rz(2.1378822) q[0];
rz(-pi) q[1];
rz(1.6414406) q[2];
sx q[2];
rz(-1.4817258) q[2];
sx q[2];
rz(-1.8235109) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2455403) q[1];
sx q[1];
rz(-2.8421092) q[1];
sx q[1];
rz(-0.36548945) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3699371) q[3];
sx q[3];
rz(-1.4691741) q[3];
sx q[3];
rz(1.9664498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6094531) q[2];
sx q[2];
rz(-0.9708465) q[2];
sx q[2];
rz(-0.54507059) q[2];
rz(0.80101454) q[3];
sx q[3];
rz(-0.89371926) q[3];
sx q[3];
rz(0.092122294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6545749) q[0];
sx q[0];
rz(-1.6681404) q[0];
sx q[0];
rz(-0.45904485) q[0];
rz(1.1162988) q[1];
sx q[1];
rz(-0.79367677) q[1];
sx q[1];
rz(-0.8265411) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7918043) q[0];
sx q[0];
rz(-1.9902181) q[0];
sx q[0];
rz(2.7546407) q[0];
x q[1];
rz(3.0439348) q[2];
sx q[2];
rz(-1.6689166) q[2];
sx q[2];
rz(-2.1922534) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.53163995) q[1];
sx q[1];
rz(-2.5594257) q[1];
sx q[1];
rz(1.9332631) q[1];
rz(1.4524121) q[3];
sx q[3];
rz(-1.9031798) q[3];
sx q[3];
rz(0.82086241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.35215846) q[2];
sx q[2];
rz(-2.9106079) q[2];
sx q[2];
rz(-2.3918772) q[2];
rz(-2.0189144) q[3];
sx q[3];
rz(-1.9176982) q[3];
sx q[3];
rz(0.96021715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66626755) q[0];
sx q[0];
rz(-0.98639494) q[0];
sx q[0];
rz(0.112003) q[0];
rz(2.9169967) q[1];
sx q[1];
rz(-1.9738395) q[1];
sx q[1];
rz(2.3602643) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3488779) q[0];
sx q[0];
rz(-1.9747866) q[0];
sx q[0];
rz(0.029026194) q[0];
rz(-0.33022837) q[2];
sx q[2];
rz(-1.4817186) q[2];
sx q[2];
rz(-1.1116127) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.427436) q[1];
sx q[1];
rz(-0.86227741) q[1];
sx q[1];
rz(-2.5664213) q[1];
rz(0.97815255) q[3];
sx q[3];
rz(-1.7053805) q[3];
sx q[3];
rz(-0.59813598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8319228) q[2];
sx q[2];
rz(-0.90193844) q[2];
sx q[2];
rz(-0.93264467) q[2];
rz(-1.0798838) q[3];
sx q[3];
rz(-2.6974758) q[3];
sx q[3];
rz(-3.1058969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.90750736) q[0];
sx q[0];
rz(-2.0103173) q[0];
sx q[0];
rz(1.750741) q[0];
rz(0.33175173) q[1];
sx q[1];
rz(-1.2650047) q[1];
sx q[1];
rz(-1.5914241) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6002613) q[0];
sx q[0];
rz(-0.92841599) q[0];
sx q[0];
rz(0.043902573) q[0];
rz(-1.5339025) q[2];
sx q[2];
rz(-1.1114745) q[2];
sx q[2];
rz(1.5806035) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0447213) q[1];
sx q[1];
rz(-1.7270178) q[1];
sx q[1];
rz(-3.0169283) q[1];
rz(-pi) q[2];
rz(0.8121536) q[3];
sx q[3];
rz(-0.5049754) q[3];
sx q[3];
rz(-0.92588378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3127689) q[2];
sx q[2];
rz(-2.2569816) q[2];
sx q[2];
rz(1.620232) q[2];
rz(-0.96794266) q[3];
sx q[3];
rz(-1.5827551) q[3];
sx q[3];
rz(-2.2255285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62436002) q[0];
sx q[0];
rz(-0.42651287) q[0];
sx q[0];
rz(-2.970001) q[0];
rz(2.102237) q[1];
sx q[1];
rz(-1.8828705) q[1];
sx q[1];
rz(-1.7003869) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0943308) q[0];
sx q[0];
rz(-0.34593317) q[0];
sx q[0];
rz(-0.37389116) q[0];
rz(-pi) q[1];
x q[1];
rz(0.19564678) q[2];
sx q[2];
rz(-0.52426978) q[2];
sx q[2];
rz(1.9877246) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4644147) q[1];
sx q[1];
rz(-0.59301335) q[1];
sx q[1];
rz(-0.89927425) q[1];
rz(1.3492382) q[3];
sx q[3];
rz(-2.5668813) q[3];
sx q[3];
rz(-1.2971085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8695996) q[2];
sx q[2];
rz(-1.7337493) q[2];
sx q[2];
rz(2.5970411) q[2];
rz(-0.49992391) q[3];
sx q[3];
rz(-1.0285503) q[3];
sx q[3];
rz(0.47206363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2917824) q[0];
sx q[0];
rz(-2.6047459) q[0];
sx q[0];
rz(0.52919069) q[0];
rz(-0.57580194) q[1];
sx q[1];
rz(-1.2628097) q[1];
sx q[1];
rz(1.1189438) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3854348) q[0];
sx q[0];
rz(-1.5498112) q[0];
sx q[0];
rz(-0.051647112) q[0];
rz(-0.37092692) q[2];
sx q[2];
rz(-1.7417522) q[2];
sx q[2];
rz(-1.8425187) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.24209472) q[1];
sx q[1];
rz(-2.3589954) q[1];
sx q[1];
rz(2.461754) q[1];
rz(1.921245) q[3];
sx q[3];
rz(-0.5548889) q[3];
sx q[3];
rz(0.12400907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.90073663) q[2];
sx q[2];
rz(-2.6764328) q[2];
sx q[2];
rz(0.71211234) q[2];
rz(-1.6093048) q[3];
sx q[3];
rz(-2.1963547) q[3];
sx q[3];
rz(1.7942662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79494548) q[0];
sx q[0];
rz(-2.0136588) q[0];
sx q[0];
rz(1.0711063) q[0];
rz(-1.2062997) q[1];
sx q[1];
rz(-1.0164398) q[1];
sx q[1];
rz(1.4869022) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7506055) q[0];
sx q[0];
rz(-1.75754) q[0];
sx q[0];
rz(2.4658527) q[0];
x q[1];
rz(-0.73545154) q[2];
sx q[2];
rz(-2.0273773) q[2];
sx q[2];
rz(-0.18010715) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7858929) q[1];
sx q[1];
rz(-1.3972056) q[1];
sx q[1];
rz(0.29533556) q[1];
rz(2.8603691) q[3];
sx q[3];
rz(-1.1373925) q[3];
sx q[3];
rz(3.0126257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6335166) q[2];
sx q[2];
rz(-0.8997007) q[2];
sx q[2];
rz(-0.077795204) q[2];
rz(2.9142694) q[3];
sx q[3];
rz(-1.1319755) q[3];
sx q[3];
rz(2.0859065) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8530497) q[0];
sx q[0];
rz(-2.396614) q[0];
sx q[0];
rz(-2.7689834) q[0];
rz(-2.695072) q[1];
sx q[1];
rz(-1.4563072) q[1];
sx q[1];
rz(2.2850697) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83607446) q[0];
sx q[0];
rz(-1.627632) q[0];
sx q[0];
rz(-1.6147805) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0356016) q[2];
sx q[2];
rz(-1.4376831) q[2];
sx q[2];
rz(3.1337381) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3831272) q[1];
sx q[1];
rz(-1.1735608) q[1];
sx q[1];
rz(0.12786156) q[1];
rz(-pi) q[2];
rz(2.4803745) q[3];
sx q[3];
rz(-2.2228129) q[3];
sx q[3];
rz(0.87424247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.38241688) q[2];
sx q[2];
rz(-1.0512115) q[2];
sx q[2];
rz(-0.55553931) q[2];
rz(0.38604745) q[3];
sx q[3];
rz(-1.6470393) q[3];
sx q[3];
rz(-2.4773795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1356708) q[0];
sx q[0];
rz(-1.6914524) q[0];
sx q[0];
rz(1.2941262) q[0];
rz(-2.338943) q[1];
sx q[1];
rz(-1.4603271) q[1];
sx q[1];
rz(-2.7535798) q[1];
rz(-0.35292179) q[2];
sx q[2];
rz(-1.6833932) q[2];
sx q[2];
rz(2.0298454) q[2];
rz(-0.84011806) q[3];
sx q[3];
rz(-1.5185322) q[3];
sx q[3];
rz(2.3867859) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
