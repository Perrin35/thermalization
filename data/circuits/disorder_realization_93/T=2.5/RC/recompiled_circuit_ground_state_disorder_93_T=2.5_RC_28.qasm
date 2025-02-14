OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.941551) q[0];
sx q[0];
rz(-0.18414944) q[0];
sx q[0];
rz(1.7848236) q[0];
rz(1.4461203) q[1];
sx q[1];
rz(-1.5674011) q[1];
sx q[1];
rz(0.27624908) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7295099) q[0];
sx q[0];
rz(-0.19664054) q[0];
sx q[0];
rz(-1.4850519) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0577781) q[2];
sx q[2];
rz(-3.028125) q[2];
sx q[2];
rz(2.6423315) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6522112) q[1];
sx q[1];
rz(-1.6141574) q[1];
sx q[1];
rz(1.5866337) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9570691) q[3];
sx q[3];
rz(-0.34466668) q[3];
sx q[3];
rz(-2.7027948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9487298) q[2];
sx q[2];
rz(-2.7728045) q[2];
sx q[2];
rz(0.31664872) q[2];
rz(-0.51490319) q[3];
sx q[3];
rz(-3.1334936) q[3];
sx q[3];
rz(0.81075877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1279476) q[0];
sx q[0];
rz(-2.7982496) q[0];
sx q[0];
rz(-3.1105594) q[0];
rz(-1.4399928) q[1];
sx q[1];
rz(-2.3411935) q[1];
sx q[1];
rz(1.4061141) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8588327) q[0];
sx q[0];
rz(-1.7525391) q[0];
sx q[0];
rz(-1.1434929) q[0];
x q[1];
rz(-3.1028264) q[2];
sx q[2];
rz(-1.3666461) q[2];
sx q[2];
rz(0.37224712) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5760562) q[1];
sx q[1];
rz(-1.4829552) q[1];
sx q[1];
rz(2.9599543) q[1];
x q[2];
rz(-0.77649103) q[3];
sx q[3];
rz(-1.463895) q[3];
sx q[3];
rz(0.85437802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7180353) q[2];
sx q[2];
rz(-0.2090629) q[2];
sx q[2];
rz(2.5162856) q[2];
rz(1.5705671) q[3];
sx q[3];
rz(-2.8795241) q[3];
sx q[3];
rz(-0.03381332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5728979) q[0];
sx q[0];
rz(-2.608572) q[0];
sx q[0];
rz(1.1385338) q[0];
rz(2.2472527) q[1];
sx q[1];
rz(-3.1411451) q[1];
sx q[1];
rz(-2.0785887) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62378685) q[0];
sx q[0];
rz(-1.8142801) q[0];
sx q[0];
rz(-2.2372449) q[0];
x q[1];
rz(-2.9350566) q[2];
sx q[2];
rz(-0.75876941) q[2];
sx q[2];
rz(-1.9643003) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0822844) q[1];
sx q[1];
rz(-1.7011146) q[1];
sx q[1];
rz(0.072029217) q[1];
x q[2];
rz(-2.390326) q[3];
sx q[3];
rz(-0.29956529) q[3];
sx q[3];
rz(-1.9825767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2669617) q[2];
sx q[2];
rz(-2.3571372) q[2];
sx q[2];
rz(-2.0060914) q[2];
rz(-1.9973495) q[3];
sx q[3];
rz(-0.091322986) q[3];
sx q[3];
rz(1.7448366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8046232) q[0];
sx q[0];
rz(-0.40358821) q[0];
sx q[0];
rz(-2.023977) q[0];
rz(-0.090967372) q[1];
sx q[1];
rz(-3.1035943) q[1];
sx q[1];
rz(1.5100347) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26474947) q[0];
sx q[0];
rz(-2.5270945) q[0];
sx q[0];
rz(-1.9351929) q[0];
rz(-pi) q[1];
rz(0.48520605) q[2];
sx q[2];
rz(-0.23113568) q[2];
sx q[2];
rz(-0.71454465) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3561894) q[1];
sx q[1];
rz(-1.265824) q[1];
sx q[1];
rz(3.0126291) q[1];
x q[2];
rz(2.1104061) q[3];
sx q[3];
rz(-2.1114719) q[3];
sx q[3];
rz(3.0669406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9522004) q[2];
sx q[2];
rz(-0.11710937) q[2];
sx q[2];
rz(-1.2516775) q[2];
rz(0.39024726) q[3];
sx q[3];
rz(-3.1249505) q[3];
sx q[3];
rz(0.19002953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8613043) q[0];
sx q[0];
rz(-2.5087293) q[0];
sx q[0];
rz(1.8508152) q[0];
rz(-3.1351008) q[1];
sx q[1];
rz(-0.26632729) q[1];
sx q[1];
rz(1.3351306) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9296918) q[0];
sx q[0];
rz(-1.9121355) q[0];
sx q[0];
rz(1.4262188) q[0];
x q[1];
rz(-0.53369227) q[2];
sx q[2];
rz(-1.6385363) q[2];
sx q[2];
rz(2.2725032) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8346429) q[1];
sx q[1];
rz(-1.6942624) q[1];
sx q[1];
rz(-0.21194063) q[1];
rz(-pi) q[2];
rz(-1.9674997) q[3];
sx q[3];
rz(-2.6894427) q[3];
sx q[3];
rz(-0.74739425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.409965) q[2];
sx q[2];
rz(-2.9755287) q[2];
sx q[2];
rz(-0.49327332) q[2];
rz(2.0491681) q[3];
sx q[3];
rz(-0.0041996669) q[3];
sx q[3];
rz(0.56001979) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6888371) q[0];
sx q[0];
rz(-1.9552564) q[0];
sx q[0];
rz(1.3415927) q[0];
rz(-1.5780559) q[1];
sx q[1];
rz(-0.23461239) q[1];
sx q[1];
rz(-2.2685307) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7811383) q[0];
sx q[0];
rz(-0.53780452) q[0];
sx q[0];
rz(-0.27307684) q[0];
x q[1];
rz(2.701917) q[2];
sx q[2];
rz(-1.5182677) q[2];
sx q[2];
rz(0.66212468) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0256589) q[1];
sx q[1];
rz(-1.6928937) q[1];
sx q[1];
rz(3.0662838) q[1];
rz(-pi) q[2];
rz(-0.21960658) q[3];
sx q[3];
rz(-1.8320932) q[3];
sx q[3];
rz(-1.1231339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0688429) q[2];
sx q[2];
rz(-1.9708956) q[2];
sx q[2];
rz(-2.0654616) q[2];
rz(0.38232803) q[3];
sx q[3];
rz(-3.101442) q[3];
sx q[3];
rz(2.373915) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1165987) q[0];
sx q[0];
rz(-0.11302639) q[0];
sx q[0];
rz(2.9207927) q[0];
rz(-2.2258672) q[1];
sx q[1];
rz(-2.974739) q[1];
sx q[1];
rz(1.487026) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0512569) q[0];
sx q[0];
rz(-1.1658586) q[0];
sx q[0];
rz(1.9693518) q[0];
x q[1];
rz(-2.0528232) q[2];
sx q[2];
rz(-2.7332167) q[2];
sx q[2];
rz(-2.0118304) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4055192) q[1];
sx q[1];
rz(-2.2637475) q[1];
sx q[1];
rz(-0.939721) q[1];
rz(3.0209157) q[3];
sx q[3];
rz(-1.4298222) q[3];
sx q[3];
rz(1.3779061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.23358341) q[2];
sx q[2];
rz(-0.073069409) q[2];
sx q[2];
rz(-2.0230587) q[2];
rz(1.3633049) q[3];
sx q[3];
rz(-0.076863591) q[3];
sx q[3];
rz(1.793687) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3240647) q[0];
sx q[0];
rz(-2.9370566) q[0];
sx q[0];
rz(-1.2943785) q[0];
rz(1.4140465) q[1];
sx q[1];
rz(-0.050948016) q[1];
sx q[1];
rz(2.9185413) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2411925) q[0];
sx q[0];
rz(-2.0471388) q[0];
sx q[0];
rz(-2.4621055) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.18977087) q[2];
sx q[2];
rz(-1.7061573) q[2];
sx q[2];
rz(3.1369097) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2168426) q[1];
sx q[1];
rz(-2.1696058) q[1];
sx q[1];
rz(-2.6393165) q[1];
x q[2];
rz(-0.87568385) q[3];
sx q[3];
rz(-1.4804262) q[3];
sx q[3];
rz(0.68708693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.53905067) q[2];
sx q[2];
rz(-1.9571303) q[2];
sx q[2];
rz(-0.43262532) q[2];
rz(-1.0299579) q[3];
sx q[3];
rz(-3.0410788) q[3];
sx q[3];
rz(-2.4242134) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3887727) q[0];
sx q[0];
rz(-2.6322375) q[0];
sx q[0];
rz(2.9260522) q[0];
rz(-2.0095677) q[1];
sx q[1];
rz(-2.3339381) q[1];
sx q[1];
rz(-1.4707123) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12080321) q[0];
sx q[0];
rz(-1.8834347) q[0];
sx q[0];
rz(-2.1759869) q[0];
x q[1];
rz(1.1342513) q[2];
sx q[2];
rz(-1.9160532) q[2];
sx q[2];
rz(2.3514812) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6214257) q[1];
sx q[1];
rz(-1.6257207) q[1];
sx q[1];
rz(-1.5406233) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.4827667) q[3];
sx q[3];
rz(-1.1454057) q[3];
sx q[3];
rz(-1.1541884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1620764) q[2];
sx q[2];
rz(-0.039204892) q[2];
sx q[2];
rz(-2.6230679) q[2];
rz(0.75820747) q[3];
sx q[3];
rz(-0.77986217) q[3];
sx q[3];
rz(-1.6583748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46490797) q[0];
sx q[0];
rz(-1.2228891) q[0];
sx q[0];
rz(1.884888) q[0];
rz(2.8769809) q[1];
sx q[1];
rz(-0.0022609641) q[1];
sx q[1];
rz(-1.9432) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6664724) q[0];
sx q[0];
rz(-1.571093) q[0];
sx q[0];
rz(-3.1350242) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9047059) q[2];
sx q[2];
rz(-1.5057202) q[2];
sx q[2];
rz(1.8462311) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8598946) q[1];
sx q[1];
rz(-1.5844189) q[1];
sx q[1];
rz(2.9071684) q[1];
rz(-pi) q[2];
rz(1.352574) q[3];
sx q[3];
rz(-0.17328158) q[3];
sx q[3];
rz(-2.1715783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4635224) q[2];
sx q[2];
rz(-3.1375725) q[2];
sx q[2];
rz(0.013997495) q[2];
rz(0.30063453) q[3];
sx q[3];
rz(-3.1379353) q[3];
sx q[3];
rz(1.5699056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8138301) q[0];
sx q[0];
rz(-1.7131085) q[0];
sx q[0];
rz(-1.499927) q[0];
rz(-3.1060524) q[1];
sx q[1];
rz(-0.4353558) q[1];
sx q[1];
rz(0.33886649) q[1];
rz(-1.8076996) q[2];
sx q[2];
rz(-2.2541004) q[2];
sx q[2];
rz(-2.83082) q[2];
rz(1.3627014) q[3];
sx q[3];
rz(-2.3632635) q[3];
sx q[3];
rz(-2.6679039) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
