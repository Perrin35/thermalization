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
rz(-2.3081245) q[0];
sx q[0];
rz(-2.3279738) q[0];
sx q[0];
rz(-1.5746434) q[0];
rz(-1.1235224) q[1];
sx q[1];
rz(-0.82668537) q[1];
sx q[1];
rz(0.5067265) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6637092) q[0];
sx q[0];
rz(-1.0179011) q[0];
sx q[0];
rz(1.5922484) q[0];
x q[1];
rz(2.0000524) q[2];
sx q[2];
rz(-0.55396336) q[2];
sx q[2];
rz(2.9316154) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6804569) q[1];
sx q[1];
rz(-1.3142921) q[1];
sx q[1];
rz(1.6272904) q[1];
rz(-2.2351609) q[3];
sx q[3];
rz(-1.6015823) q[3];
sx q[3];
rz(-0.258729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7099521) q[2];
sx q[2];
rz(-1.9441354) q[2];
sx q[2];
rz(0.54090995) q[2];
rz(1.3276118) q[3];
sx q[3];
rz(-1.6874467) q[3];
sx q[3];
rz(2.4732164) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3427908) q[0];
sx q[0];
rz(-1.3725766) q[0];
sx q[0];
rz(-2.0846682) q[0];
rz(-2.7150555) q[1];
sx q[1];
rz(-1.6273472) q[1];
sx q[1];
rz(0.0043409745) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98470178) q[0];
sx q[0];
rz(-0.92388573) q[0];
sx q[0];
rz(2.6527651) q[0];
rz(2.1059733) q[2];
sx q[2];
rz(-0.7053203) q[2];
sx q[2];
rz(-2.3065604) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7340034) q[1];
sx q[1];
rz(-0.079691039) q[1];
sx q[1];
rz(1.4400287) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6326174) q[3];
sx q[3];
rz(-1.0366716) q[3];
sx q[3];
rz(1.5388863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1987552) q[2];
sx q[2];
rz(-1.402907) q[2];
sx q[2];
rz(2.7665561) q[2];
rz(0.047920553) q[3];
sx q[3];
rz(-1.5791357) q[3];
sx q[3];
rz(-2.9722049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7584943) q[0];
sx q[0];
rz(-2.5653745) q[0];
sx q[0];
rz(-1.5118442) q[0];
rz(2.0747298) q[1];
sx q[1];
rz(-1.8424415) q[1];
sx q[1];
rz(-2.3323434) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3602995) q[0];
sx q[0];
rz(-1.5341833) q[0];
sx q[0];
rz(3.1017324) q[0];
rz(-pi) q[1];
rz(0.079945076) q[2];
sx q[2];
rz(-2.4221651) q[2];
sx q[2];
rz(-1.5647174) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7429333) q[1];
sx q[1];
rz(-2.3509532) q[1];
sx q[1];
rz(-1.8439951) q[1];
x q[2];
rz(-2.7797805) q[3];
sx q[3];
rz(-1.5734976) q[3];
sx q[3];
rz(2.5914374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2201436) q[2];
sx q[2];
rz(-2.4704832) q[2];
sx q[2];
rz(-1.6371833) q[2];
rz(2.1680221) q[3];
sx q[3];
rz(-0.74876553) q[3];
sx q[3];
rz(2.0481295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.022738986) q[0];
sx q[0];
rz(-1.5819419) q[0];
sx q[0];
rz(2.7795025) q[0];
rz(1.7936961) q[1];
sx q[1];
rz(-1.3840414) q[1];
sx q[1];
rz(-1.0628343) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75116144) q[0];
sx q[0];
rz(-1.4289723) q[0];
sx q[0];
rz(-1.7447937) q[0];
rz(-pi) q[1];
rz(0.21740977) q[2];
sx q[2];
rz(-0.98071456) q[2];
sx q[2];
rz(-1.3979531) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.98403) q[1];
sx q[1];
rz(-1.8631991) q[1];
sx q[1];
rz(-1.2127905) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.70433299) q[3];
sx q[3];
rz(-2.8129431) q[3];
sx q[3];
rz(2.8002594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7245543) q[2];
sx q[2];
rz(-0.78196708) q[2];
sx q[2];
rz(-0.057272591) q[2];
rz(-0.96632424) q[3];
sx q[3];
rz(-0.94760197) q[3];
sx q[3];
rz(-1.2053325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9178977) q[0];
sx q[0];
rz(-2.2530093) q[0];
sx q[0];
rz(-2.5433871) q[0];
rz(-1.5679476) q[1];
sx q[1];
rz(-1.4718098) q[1];
sx q[1];
rz(-2.9817458) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.689834) q[0];
sx q[0];
rz(-1.7421113) q[0];
sx q[0];
rz(2.3687045) q[0];
rz(-pi) q[1];
rz(2.9303126) q[2];
sx q[2];
rz(-1.5108372) q[2];
sx q[2];
rz(0.8983537) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5626403) q[1];
sx q[1];
rz(-0.75030164) q[1];
sx q[1];
rz(-0.57344389) q[1];
rz(-1.615404) q[3];
sx q[3];
rz(-1.6470223) q[3];
sx q[3];
rz(-1.7895229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9827031) q[2];
sx q[2];
rz(-2.7003459) q[2];
sx q[2];
rz(-0.59054792) q[2];
rz(-2.1741137) q[3];
sx q[3];
rz(-1.5524105) q[3];
sx q[3];
rz(-1.9578741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9443053) q[0];
sx q[0];
rz(-0.66216457) q[0];
sx q[0];
rz(-1.039132) q[0];
rz(-2.2100673) q[1];
sx q[1];
rz(-0.68382278) q[1];
sx q[1];
rz(-1.1087803) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1241348) q[0];
sx q[0];
rz(-2.0724247) q[0];
sx q[0];
rz(-2.5162016) q[0];
x q[1];
rz(-0.81056912) q[2];
sx q[2];
rz(-2.2165073) q[2];
sx q[2];
rz(-3.1217074) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.30516827) q[1];
sx q[1];
rz(-2.1427611) q[1];
sx q[1];
rz(1.9981153) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.118947) q[3];
sx q[3];
rz(-1.7860869) q[3];
sx q[3];
rz(-1.5529274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.335677) q[2];
sx q[2];
rz(-0.47199619) q[2];
sx q[2];
rz(2.8373888) q[2];
rz(0.16600569) q[3];
sx q[3];
rz(-1.3622354) q[3];
sx q[3];
rz(-1.5267052) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6840927) q[0];
sx q[0];
rz(-0.2245716) q[0];
sx q[0];
rz(-1.7042879) q[0];
rz(2.8999088) q[1];
sx q[1];
rz(-1.4983404) q[1];
sx q[1];
rz(2.2041352) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3110549) q[0];
sx q[0];
rz(-1.8094581) q[0];
sx q[0];
rz(2.412459) q[0];
rz(-1.2901914) q[2];
sx q[2];
rz(-1.5414503) q[2];
sx q[2];
rz(1.2236694) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7852206) q[1];
sx q[1];
rz(-0.8360148) q[1];
sx q[1];
rz(-2.4874175) q[1];
rz(2.695153) q[3];
sx q[3];
rz(-1.0390128) q[3];
sx q[3];
rz(1.4908294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.23718701) q[2];
sx q[2];
rz(-2.0439456) q[2];
sx q[2];
rz(-0.22906765) q[2];
rz(-1.9704874) q[3];
sx q[3];
rz(-1.7541211) q[3];
sx q[3];
rz(-2.4598725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1218629) q[0];
sx q[0];
rz(-0.64058146) q[0];
sx q[0];
rz(-1.4135452) q[0];
rz(1.9238663) q[1];
sx q[1];
rz(-1.7037337) q[1];
sx q[1];
rz(-2.6103643) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41252688) q[0];
sx q[0];
rz(-0.382889) q[0];
sx q[0];
rz(-1.0273496) q[0];
rz(-2.7338303) q[2];
sx q[2];
rz(-0.86155781) q[2];
sx q[2];
rz(0.10283486) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5672877) q[1];
sx q[1];
rz(-2.0101846) q[1];
sx q[1];
rz(-0.69171016) q[1];
x q[2];
rz(-1.2770428) q[3];
sx q[3];
rz(-1.2790907) q[3];
sx q[3];
rz(-0.13495378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.56500498) q[2];
sx q[2];
rz(-1.7030623) q[2];
sx q[2];
rz(-2.3459404) q[2];
rz(-2.3203881) q[3];
sx q[3];
rz(-0.20902769) q[3];
sx q[3];
rz(-2.8388099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8177556) q[0];
sx q[0];
rz(-0.19446401) q[0];
sx q[0];
rz(1.6545464) q[0];
rz(-1.6200292) q[1];
sx q[1];
rz(-1.227102) q[1];
sx q[1];
rz(-1.0618658) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8811765) q[0];
sx q[0];
rz(-1.9435391) q[0];
sx q[0];
rz(2.580852) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8650811) q[2];
sx q[2];
rz(-1.6194034) q[2];
sx q[2];
rz(-0.44845861) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3652851) q[1];
sx q[1];
rz(-2.690491) q[1];
sx q[1];
rz(1.9264313) q[1];
rz(-1.1561772) q[3];
sx q[3];
rz(-1.7545757) q[3];
sx q[3];
rz(-1.444425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1294641) q[2];
sx q[2];
rz(-2.0120554) q[2];
sx q[2];
rz(3.0322292) q[2];
rz(-2.8520285) q[3];
sx q[3];
rz(-1.7508296) q[3];
sx q[3];
rz(-2.4880828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7107723) q[0];
sx q[0];
rz(-1.9662974) q[0];
sx q[0];
rz(1.1547422) q[0];
rz(0.42354241) q[1];
sx q[1];
rz(-1.4332899) q[1];
sx q[1];
rz(0.67798859) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38289878) q[0];
sx q[0];
rz(-2.3385171) q[0];
sx q[0];
rz(2.4944958) q[0];
rz(-1.7977436) q[2];
sx q[2];
rz(-1.316202) q[2];
sx q[2];
rz(0.4190906) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0748234) q[1];
sx q[1];
rz(-0.69157234) q[1];
sx q[1];
rz(-2.2645386) q[1];
rz(-pi) q[2];
rz(0.20991169) q[3];
sx q[3];
rz(-0.60214199) q[3];
sx q[3];
rz(-1.9937552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.28422156) q[2];
sx q[2];
rz(-0.81615964) q[2];
sx q[2];
rz(-0.59874272) q[2];
rz(1.6019999) q[3];
sx q[3];
rz(-1.2076104) q[3];
sx q[3];
rz(-2.1046751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4620062) q[0];
sx q[0];
rz(-2.7988667) q[0];
sx q[0];
rz(-1.3890247) q[0];
rz(3.1287843) q[1];
sx q[1];
rz(-0.3777596) q[1];
sx q[1];
rz(-1.6526745) q[1];
rz(2.2118582) q[2];
sx q[2];
rz(-1.6597181) q[2];
sx q[2];
rz(1.5372288) q[2];
rz(-1.5725224) q[3];
sx q[3];
rz(-0.55894096) q[3];
sx q[3];
rz(-0.10645549) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
