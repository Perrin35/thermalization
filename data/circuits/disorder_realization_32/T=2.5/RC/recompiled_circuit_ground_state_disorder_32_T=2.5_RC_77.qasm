OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.6736421) q[0];
sx q[0];
rz(7.5543348) q[0];
sx q[0];
rz(13.332097) q[0];
rz(-0.39628059) q[1];
sx q[1];
rz(-0.09274617) q[1];
sx q[1];
rz(0.5527817) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6792949) q[0];
sx q[0];
rz(-2.5640985) q[0];
sx q[0];
rz(-0.80120835) q[0];
rz(2.4713995) q[2];
sx q[2];
rz(-2.2171221) q[2];
sx q[2];
rz(2.0984416) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1172792) q[1];
sx q[1];
rz(-1.2543884) q[1];
sx q[1];
rz(-0.98543075) q[1];
x q[2];
rz(0.90867282) q[3];
sx q[3];
rz(-1.931802) q[3];
sx q[3];
rz(2.8761325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7736241) q[2];
sx q[2];
rz(-1.6398733) q[2];
sx q[2];
rz(3.0527414) q[2];
rz(0.39696524) q[3];
sx q[3];
rz(-1.0395972) q[3];
sx q[3];
rz(1.6840434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48001584) q[0];
sx q[0];
rz(-0.50742298) q[0];
sx q[0];
rz(0.85451025) q[0];
rz(0.62141934) q[1];
sx q[1];
rz(-2.8050551) q[1];
sx q[1];
rz(-2.0889166) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54224211) q[0];
sx q[0];
rz(-0.17048888) q[0];
sx q[0];
rz(0.64298274) q[0];
rz(0.82727716) q[2];
sx q[2];
rz(-0.57663871) q[2];
sx q[2];
rz(-1.0232384) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4667644) q[1];
sx q[1];
rz(-0.75777868) q[1];
sx q[1];
rz(-1.8653052) q[1];
rz(1.5684116) q[3];
sx q[3];
rz(-0.60394507) q[3];
sx q[3];
rz(1.2177474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4738327) q[2];
sx q[2];
rz(-0.86060539) q[2];
sx q[2];
rz(0.94998002) q[2];
rz(1.0085603) q[3];
sx q[3];
rz(-0.63298321) q[3];
sx q[3];
rz(-0.27957255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.1968593) q[0];
sx q[0];
rz(-1.9254528) q[0];
sx q[0];
rz(1.3470294) q[0];
rz(2.0836209) q[1];
sx q[1];
rz(-2.8949013) q[1];
sx q[1];
rz(3.0314441) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63907097) q[0];
sx q[0];
rz(-2.5701231) q[0];
sx q[0];
rz(-2.1854464) q[0];
rz(1.4227432) q[2];
sx q[2];
rz(-1.4880991) q[2];
sx q[2];
rz(-2.952864) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.55787841) q[1];
sx q[1];
rz(-1.4643351) q[1];
sx q[1];
rz(2.540968) q[1];
rz(-pi) q[2];
rz(-2.0581117) q[3];
sx q[3];
rz(-1.1705453) q[3];
sx q[3];
rz(2.6292173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9614253) q[2];
sx q[2];
rz(-1.578293) q[2];
sx q[2];
rz(0.041672826) q[2];
rz(2.9336119) q[3];
sx q[3];
rz(-2.9198923) q[3];
sx q[3];
rz(-0.25377932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5792849) q[0];
sx q[0];
rz(-1.8203745) q[0];
sx q[0];
rz(-0.98096171) q[0];
rz(-2.7085069) q[1];
sx q[1];
rz(-1.4739477) q[1];
sx q[1];
rz(-1.6281737) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53674066) q[0];
sx q[0];
rz(-1.3910595) q[0];
sx q[0];
rz(-2.7374703) q[0];
x q[1];
rz(3.1279706) q[2];
sx q[2];
rz(-2.2421466) q[2];
sx q[2];
rz(1.4289795) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.42964307) q[1];
sx q[1];
rz(-2.6338661) q[1];
sx q[1];
rz(1.7832548) q[1];
rz(-pi) q[2];
rz(-1.9497197) q[3];
sx q[3];
rz(-2.4746568) q[3];
sx q[3];
rz(-0.88465026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.607434) q[2];
sx q[2];
rz(-0.51924339) q[2];
sx q[2];
rz(-2.329211) q[2];
rz(1.8937998) q[3];
sx q[3];
rz(-1.8915853) q[3];
sx q[3];
rz(-1.2468503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(1.8060551) q[0];
sx q[0];
rz(-2.1195109) q[0];
sx q[0];
rz(1.6988276) q[0];
rz(0.45817786) q[1];
sx q[1];
rz(-1.5135601) q[1];
sx q[1];
rz(2.3853669) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0201705) q[0];
sx q[0];
rz(-2.4947632) q[0];
sx q[0];
rz(1.8298762) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.67839038) q[2];
sx q[2];
rz(-0.76324651) q[2];
sx q[2];
rz(2.2445172) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1444662) q[1];
sx q[1];
rz(-1.6580075) q[1];
sx q[1];
rz(3.0060796) q[1];
rz(-pi) q[2];
rz(1.9147768) q[3];
sx q[3];
rz(-1.5369253) q[3];
sx q[3];
rz(-0.86465166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4592287) q[2];
sx q[2];
rz(-1.3882779) q[2];
sx q[2];
rz(0.96561042) q[2];
rz(-2.5718578) q[3];
sx q[3];
rz(-2.0163048) q[3];
sx q[3];
rz(-1.4558571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4023034) q[0];
sx q[0];
rz(-1.1981413) q[0];
sx q[0];
rz(-1.7560316) q[0];
rz(-0.7631453) q[1];
sx q[1];
rz(-1.6940073) q[1];
sx q[1];
rz(-0.10173434) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4674326) q[0];
sx q[0];
rz(-1.4786353) q[0];
sx q[0];
rz(-2.1919495) q[0];
rz(2.7587682) q[2];
sx q[2];
rz(-1.6144132) q[2];
sx q[2];
rz(-0.62784615) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2925955) q[1];
sx q[1];
rz(-1.6319872) q[1];
sx q[1];
rz(-1.2996177) q[1];
rz(-pi) q[2];
rz(-0.88512986) q[3];
sx q[3];
rz(-0.99381522) q[3];
sx q[3];
rz(0.5761522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0330641) q[2];
sx q[2];
rz(-1.6988924) q[2];
sx q[2];
rz(0.70890439) q[2];
rz(2.3809643) q[3];
sx q[3];
rz(-1.9899188) q[3];
sx q[3];
rz(-2.2061548) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0660504) q[0];
sx q[0];
rz(-1.9882555) q[0];
sx q[0];
rz(0.72108889) q[0];
rz(-0.9779633) q[1];
sx q[1];
rz(-0.55897346) q[1];
sx q[1];
rz(1.5616547) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2943731) q[0];
sx q[0];
rz(-0.10589639) q[0];
sx q[0];
rz(0.76283331) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8992813) q[2];
sx q[2];
rz(-2.152395) q[2];
sx q[2];
rz(-2.966449) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9813919) q[1];
sx q[1];
rz(-1.0202279) q[1];
sx q[1];
rz(2.0345313) q[1];
x q[2];
rz(0.35646371) q[3];
sx q[3];
rz(-0.93178669) q[3];
sx q[3];
rz(-1.654341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8568023) q[2];
sx q[2];
rz(-0.52758104) q[2];
sx q[2];
rz(-1.6214726) q[2];
rz(2.704845) q[3];
sx q[3];
rz(-0.89478651) q[3];
sx q[3];
rz(-2.6020218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.0590416) q[0];
sx q[0];
rz(-0.83913791) q[0];
sx q[0];
rz(0.82373291) q[0];
rz(-1.1388904) q[1];
sx q[1];
rz(-1.3528119) q[1];
sx q[1];
rz(-1.9653856) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.091789548) q[0];
sx q[0];
rz(-0.50423008) q[0];
sx q[0];
rz(-1.6635899) q[0];
rz(-pi) q[1];
rz(0.55032879) q[2];
sx q[2];
rz(-2.4541353) q[2];
sx q[2];
rz(-2.2000809) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2580793) q[1];
sx q[1];
rz(-0.23830676) q[1];
sx q[1];
rz(-3.1001453) q[1];
x q[2];
rz(-3.0823364) q[3];
sx q[3];
rz(-2.5496581) q[3];
sx q[3];
rz(1.8979285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.99522432) q[2];
sx q[2];
rz(-2.5547042) q[2];
sx q[2];
rz(2.5541019) q[2];
rz(-2.7557709) q[3];
sx q[3];
rz(-2.2930175) q[3];
sx q[3];
rz(0.90252701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0226444) q[0];
sx q[0];
rz(-0.97624874) q[0];
sx q[0];
rz(2.5318085) q[0];
rz(-1.7976286) q[1];
sx q[1];
rz(-1.983843) q[1];
sx q[1];
rz(0.15300289) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5822844) q[0];
sx q[0];
rz(-1.1049718) q[0];
sx q[0];
rz(1.5927754) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.010092334) q[2];
sx q[2];
rz(-2.3463425) q[2];
sx q[2];
rz(-2.681596) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.24125946) q[1];
sx q[1];
rz(-1.5521084) q[1];
sx q[1];
rz(-2.2128723) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1243662) q[3];
sx q[3];
rz(-2.3690936) q[3];
sx q[3];
rz(1.8541418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1028221) q[2];
sx q[2];
rz(-1.5800579) q[2];
sx q[2];
rz(0.62225303) q[2];
rz(-1.2004987) q[3];
sx q[3];
rz(-1.4193204) q[3];
sx q[3];
rz(2.5672369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1822561) q[0];
sx q[0];
rz(-3.0756364) q[0];
sx q[0];
rz(0.35520735) q[0];
rz(-1.1025053) q[1];
sx q[1];
rz(-3.1246154) q[1];
sx q[1];
rz(-2.4818518) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3400354) q[0];
sx q[0];
rz(-1.7424955) q[0];
sx q[0];
rz(2.4853909) q[0];
rz(-pi) q[1];
rz(-3.0142886) q[2];
sx q[2];
rz(-1.9812756) q[2];
sx q[2];
rz(-2.8212027) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.5774571) q[1];
sx q[1];
rz(-0.45960765) q[1];
sx q[1];
rz(2.3338302) q[1];
rz(1.6361314) q[3];
sx q[3];
rz(-2.4628277) q[3];
sx q[3];
rz(2.9836054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3416662) q[2];
sx q[2];
rz(-0.081192668) q[2];
sx q[2];
rz(2.9426835) q[2];
rz(2.343446) q[3];
sx q[3];
rz(-1.9559559) q[3];
sx q[3];
rz(1.859349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9195008) q[0];
sx q[0];
rz(-1.5380479) q[0];
sx q[0];
rz(2.0684239) q[0];
rz(-0.78492289) q[1];
sx q[1];
rz(-2.2833318) q[1];
sx q[1];
rz(-2.2699184) q[1];
rz(-1.6941407) q[2];
sx q[2];
rz(-1.5236749) q[2];
sx q[2];
rz(2.9070367) q[2];
rz(-2.8002515) q[3];
sx q[3];
rz(-2.4433856) q[3];
sx q[3];
rz(-2.7563765) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
