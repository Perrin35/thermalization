OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.61092678) q[0];
sx q[0];
rz(-0.93902421) q[0];
sx q[0];
rz(-0.0052069081) q[0];
rz(-1.3448673) q[1];
sx q[1];
rz(-1.1562647) q[1];
sx q[1];
rz(-1.9519238) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0931041) q[0];
sx q[0];
rz(-0.63360533) q[0];
sx q[0];
rz(-2.5392169) q[0];
rz(0.18007937) q[2];
sx q[2];
rz(-1.4305978) q[2];
sx q[2];
rz(2.7498498) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7835044) q[1];
sx q[1];
rz(-1.7184098) q[1];
sx q[1];
rz(-1.258177) q[1];
rz(-pi) q[2];
rz(2.1055929) q[3];
sx q[3];
rz(-0.70123226) q[3];
sx q[3];
rz(2.0555156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.71620119) q[2];
sx q[2];
rz(-2.3925245) q[2];
sx q[2];
rz(2.544196) q[2];
rz(-1.776009) q[3];
sx q[3];
rz(-1.7979012) q[3];
sx q[3];
rz(-1.8610154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4269203) q[0];
sx q[0];
rz(-0.5517813) q[0];
sx q[0];
rz(2.8080217) q[0];
rz(1.0936273) q[1];
sx q[1];
rz(-0.90197864) q[1];
sx q[1];
rz(3.0283668) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9202168) q[0];
sx q[0];
rz(-2.7144103) q[0];
sx q[0];
rz(-1.6643874) q[0];
x q[1];
rz(-0.073268854) q[2];
sx q[2];
rz(-1.5904625) q[2];
sx q[2];
rz(-2.427223) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7659521) q[1];
sx q[1];
rz(-1.1728219) q[1];
sx q[1];
rz(-2.6938733) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1749288) q[3];
sx q[3];
rz(-1.767744) q[3];
sx q[3];
rz(0.40991022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1229646) q[2];
sx q[2];
rz(-2.661442) q[2];
sx q[2];
rz(-0.22228995) q[2];
rz(-2.9120581) q[3];
sx q[3];
rz(-0.73733202) q[3];
sx q[3];
rz(-0.078331746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6753321) q[0];
sx q[0];
rz(-1.5706797) q[0];
sx q[0];
rz(-2.2608742) q[0];
rz(-1.0473898) q[1];
sx q[1];
rz(-1.9347582) q[1];
sx q[1];
rz(0.12761322) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.815044) q[0];
sx q[0];
rz(-2.43988) q[0];
sx q[0];
rz(-1.2998507) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3348546) q[2];
sx q[2];
rz(-2.2819937) q[2];
sx q[2];
rz(-2.9850609) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.322791) q[1];
sx q[1];
rz(-1.6530419) q[1];
sx q[1];
rz(1.8104042) q[1];
rz(-pi) q[2];
rz(1.2283986) q[3];
sx q[3];
rz(-0.6558334) q[3];
sx q[3];
rz(-2.0369903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7814653) q[2];
sx q[2];
rz(-1.1294304) q[2];
sx q[2];
rz(-0.310251) q[2];
rz(0.34960738) q[3];
sx q[3];
rz(-0.6844371) q[3];
sx q[3];
rz(-1.413697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.67868245) q[0];
sx q[0];
rz(-1.3218198) q[0];
sx q[0];
rz(0.15790766) q[0];
rz(-0.36610106) q[1];
sx q[1];
rz(-0.78916517) q[1];
sx q[1];
rz(-2.7691832) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87293738) q[0];
sx q[0];
rz(-2.4252709) q[0];
sx q[0];
rz(1.3852081) q[0];
rz(2.1521058) q[2];
sx q[2];
rz(-1.8539691) q[2];
sx q[2];
rz(0.95209939) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6652674) q[1];
sx q[1];
rz(-1.421126) q[1];
sx q[1];
rz(-1.1900351) q[1];
rz(2.1129205) q[3];
sx q[3];
rz(-0.91915932) q[3];
sx q[3];
rz(1.2937014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7238414) q[2];
sx q[2];
rz(-2.4390742) q[2];
sx q[2];
rz(2.8692029) q[2];
rz(0.90879905) q[3];
sx q[3];
rz(-2.2464928) q[3];
sx q[3];
rz(-0.4549543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2588147) q[0];
sx q[0];
rz(-0.60084501) q[0];
sx q[0];
rz(-1.1313261) q[0];
rz(-2.5436026) q[1];
sx q[1];
rz(-2.2441041) q[1];
sx q[1];
rz(-0.67684832) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2535431) q[0];
sx q[0];
rz(-1.4918259) q[0];
sx q[0];
rz(2.700564) q[0];
rz(-2.9025181) q[2];
sx q[2];
rz(-0.74547807) q[2];
sx q[2];
rz(2.4649232) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0332898) q[1];
sx q[1];
rz(-1.3160333) q[1];
sx q[1];
rz(1.5441896) q[1];
rz(-1.5464209) q[3];
sx q[3];
rz(-1.1773603) q[3];
sx q[3];
rz(-2.3778621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.1375492) q[2];
sx q[2];
rz(-1.153959) q[2];
sx q[2];
rz(-1.1452902) q[2];
rz(2.9925313) q[3];
sx q[3];
rz(-0.59864932) q[3];
sx q[3];
rz(2.1242274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0155708) q[0];
sx q[0];
rz(-0.64783043) q[0];
sx q[0];
rz(1.6824678) q[0];
rz(0.74516621) q[1];
sx q[1];
rz(-0.35062733) q[1];
sx q[1];
rz(-0.93313342) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81290302) q[0];
sx q[0];
rz(-1.110154) q[0];
sx q[0];
rz(-2.8115389) q[0];
rz(-pi) q[1];
rz(-0.60284166) q[2];
sx q[2];
rz(-1.8800929) q[2];
sx q[2];
rz(-1.7939292) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9698525) q[1];
sx q[1];
rz(-1.9386567) q[1];
sx q[1];
rz(0.0024585558) q[1];
rz(-pi) q[2];
rz(-0.91306367) q[3];
sx q[3];
rz(-1.7115895) q[3];
sx q[3];
rz(-1.4482244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.78559819) q[2];
sx q[2];
rz(-0.22452536) q[2];
sx q[2];
rz(0.28665001) q[2];
rz(-2.8921228) q[3];
sx q[3];
rz(-1.1688066) q[3];
sx q[3];
rz(2.6732895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.048112415) q[0];
sx q[0];
rz(-1.965964) q[0];
sx q[0];
rz(1.6749143) q[0];
rz(-2.1215227) q[1];
sx q[1];
rz(-1.6721098) q[1];
sx q[1];
rz(-1.0010304) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6769343) q[0];
sx q[0];
rz(-1.5334198) q[0];
sx q[0];
rz(-1.8514762) q[0];
x q[1];
rz(1.2008576) q[2];
sx q[2];
rz(-1.2518034) q[2];
sx q[2];
rz(1.5202886) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.475358) q[1];
sx q[1];
rz(-1.959414) q[1];
sx q[1];
rz(2.9734441) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4909199) q[3];
sx q[3];
rz(-1.4015897) q[3];
sx q[3];
rz(-2.8229439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.004185685) q[2];
sx q[2];
rz(-1.9150532) q[2];
sx q[2];
rz(0.10350791) q[2];
rz(-0.59182709) q[3];
sx q[3];
rz(-1.3261565) q[3];
sx q[3];
rz(2.3039968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8822534) q[0];
sx q[0];
rz(-2.1253724) q[0];
sx q[0];
rz(-2.5296339) q[0];
rz(1.7991964) q[1];
sx q[1];
rz(-1.1609158) q[1];
sx q[1];
rz(0.38988316) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.10605) q[0];
sx q[0];
rz(-2.4379726) q[0];
sx q[0];
rz(-2.5662867) q[0];
rz(-pi) q[1];
rz(-1.9837603) q[2];
sx q[2];
rz(-1.0997084) q[2];
sx q[2];
rz(0.72380356) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0695614) q[1];
sx q[1];
rz(-1.2417292) q[1];
sx q[1];
rz(-2.0272994) q[1];
x q[2];
rz(-0.26384683) q[3];
sx q[3];
rz(-0.90875328) q[3];
sx q[3];
rz(2.2784233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8395681) q[2];
sx q[2];
rz(-2.6111157) q[2];
sx q[2];
rz(-2.4273382) q[2];
rz(-2.2438625) q[3];
sx q[3];
rz(-1.1313181) q[3];
sx q[3];
rz(2.5695661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4488688) q[0];
sx q[0];
rz(-2.7063997) q[0];
sx q[0];
rz(2.3642448) q[0];
rz(-2.2946987) q[1];
sx q[1];
rz(-2.6234026) q[1];
sx q[1];
rz(2.8651967) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9696635) q[0];
sx q[0];
rz(-1.164325) q[0];
sx q[0];
rz(-1.3967692) q[0];
rz(-pi) q[1];
x q[1];
rz(0.22819613) q[2];
sx q[2];
rz(-1.3820634) q[2];
sx q[2];
rz(-0.91918321) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8564057) q[1];
sx q[1];
rz(-1.7289843) q[1];
sx q[1];
rz(1.7941979) q[1];
x q[2];
rz(-1.8180088) q[3];
sx q[3];
rz(-2.1650378) q[3];
sx q[3];
rz(-2.0127399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2892264) q[2];
sx q[2];
rz(-1.6324531) q[2];
sx q[2];
rz(-3.0140871) q[2];
rz(-3.1048807) q[3];
sx q[3];
rz(-2.3214985) q[3];
sx q[3];
rz(-1.9071074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.4002832) q[0];
sx q[0];
rz(-2.648634) q[0];
sx q[0];
rz(2.7888443) q[0];
rz(-2.5648975) q[1];
sx q[1];
rz(-0.99027094) q[1];
sx q[1];
rz(-1.030285) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0477662) q[0];
sx q[0];
rz(-0.86137912) q[0];
sx q[0];
rz(1.8295656) q[0];
rz(-pi) q[1];
x q[1];
rz(2.455591) q[2];
sx q[2];
rz(-1.8324319) q[2];
sx q[2];
rz(0.54856448) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.87957571) q[1];
sx q[1];
rz(-2.587662) q[1];
sx q[1];
rz(-2.5074156) q[1];
rz(-1.9029721) q[3];
sx q[3];
rz(-2.1579451) q[3];
sx q[3];
rz(1.7410994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0742005) q[2];
sx q[2];
rz(-1.3211162) q[2];
sx q[2];
rz(0.34118787) q[2];
rz(-2.4009005) q[3];
sx q[3];
rz(-2.1518028) q[3];
sx q[3];
rz(-3.0497131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0108903) q[0];
sx q[0];
rz(-1.9938835) q[0];
sx q[0];
rz(2.4947517) q[0];
rz(2.4717992) q[1];
sx q[1];
rz(-2.5302946) q[1];
sx q[1];
rz(-1.5763462) q[1];
rz(-1.7651991) q[2];
sx q[2];
rz(-1.1444848) q[2];
sx q[2];
rz(-0.3111006) q[2];
rz(0.66486799) q[3];
sx q[3];
rz(-2.3330199) q[3];
sx q[3];
rz(-1.4521269) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
