OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.527737) q[0];
sx q[0];
rz(-1.4976488) q[0];
sx q[0];
rz(-2.3117476) q[0];
rz(3.9217477) q[1];
sx q[1];
rz(5.2182066) q[1];
sx q[1];
rz(10.301104) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0759461) q[0];
sx q[0];
rz(-1.3407205) q[0];
sx q[0];
rz(-2.7394657) q[0];
rz(-pi) q[1];
rz(1.285577) q[2];
sx q[2];
rz(-0.60456317) q[2];
sx q[2];
rz(2.4100458) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6688924) q[1];
sx q[1];
rz(-1.9007705) q[1];
sx q[1];
rz(-0.015923576) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.34378864) q[3];
sx q[3];
rz(-0.4292092) q[3];
sx q[3];
rz(-0.49376282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.17065389) q[2];
sx q[2];
rz(-1.8654414) q[2];
sx q[2];
rz(-0.7286287) q[2];
rz(-0.5209926) q[3];
sx q[3];
rz(-2.1803768) q[3];
sx q[3];
rz(2.9339824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.306863) q[0];
sx q[0];
rz(-1.1704209) q[0];
sx q[0];
rz(-2.0200502) q[0];
rz(-0.25575486) q[1];
sx q[1];
rz(-1.6633727) q[1];
sx q[1];
rz(2.2671525) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6073608) q[0];
sx q[0];
rz(-0.53593862) q[0];
sx q[0];
rz(2.1952654) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8057683) q[2];
sx q[2];
rz(-2.5532696) q[2];
sx q[2];
rz(2.7749643) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2494617) q[1];
sx q[1];
rz(-0.79626894) q[1];
sx q[1];
rz(-2.4888121) q[1];
x q[2];
rz(0.053697649) q[3];
sx q[3];
rz(-1.8397545) q[3];
sx q[3];
rz(2.4503436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.740739) q[2];
sx q[2];
rz(-2.6541371) q[2];
sx q[2];
rz(-0.43593105) q[2];
rz(-0.68108264) q[3];
sx q[3];
rz(-2.3705132) q[3];
sx q[3];
rz(2.7387103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9044559) q[0];
sx q[0];
rz(-2.853892) q[0];
sx q[0];
rz(-2.0667734) q[0];
rz(2.3020321) q[1];
sx q[1];
rz(-0.81937516) q[1];
sx q[1];
rz(-0.39594617) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.031035) q[0];
sx q[0];
rz(-1.1764515) q[0];
sx q[0];
rz(-2.1179384) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8648519) q[2];
sx q[2];
rz(-1.5573504) q[2];
sx q[2];
rz(-0.80317342) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.46519854) q[1];
sx q[1];
rz(-0.68323831) q[1];
sx q[1];
rz(-0.015124358) q[1];
rz(-pi) q[2];
rz(-1.5941761) q[3];
sx q[3];
rz(-1.3153207) q[3];
sx q[3];
rz(1.3083096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5376771) q[2];
sx q[2];
rz(-1.2664412) q[2];
sx q[2];
rz(-0.23920693) q[2];
rz(-0.075332969) q[3];
sx q[3];
rz(-2.0276666) q[3];
sx q[3];
rz(1.8384365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4199715) q[0];
sx q[0];
rz(-2.2890838) q[0];
sx q[0];
rz(-2.3216632) q[0];
rz(-2.6539102) q[1];
sx q[1];
rz(-2.2380424) q[1];
sx q[1];
rz(0.23342361) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6596286) q[0];
sx q[0];
rz(-0.11419645) q[0];
sx q[0];
rz(1.0114848) q[0];
x q[1];
rz(2.8308224) q[2];
sx q[2];
rz(-0.66590532) q[2];
sx q[2];
rz(-3.0096465) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9277716) q[1];
sx q[1];
rz(-1.2328316) q[1];
sx q[1];
rz(2.4073699) q[1];
x q[2];
rz(-0.78419533) q[3];
sx q[3];
rz(-1.8521706) q[3];
sx q[3];
rz(-2.1649862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0908115) q[2];
sx q[2];
rz(-1.8596785) q[2];
sx q[2];
rz(-0.74742571) q[2];
rz(2.9181972) q[3];
sx q[3];
rz(-2.5441393) q[3];
sx q[3];
rz(2.8994765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6915879) q[0];
sx q[0];
rz(-1.7084028) q[0];
sx q[0];
rz(3.0773556) q[0];
rz(-2.1977987) q[1];
sx q[1];
rz(-0.72729021) q[1];
sx q[1];
rz(-0.76104004) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9984263) q[0];
sx q[0];
rz(-0.26253065) q[0];
sx q[0];
rz(-1.4545928) q[0];
rz(-1.3752851) q[2];
sx q[2];
rz(-1.3720781) q[2];
sx q[2];
rz(2.9544601) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6619751) q[1];
sx q[1];
rz(-1.8912589) q[1];
sx q[1];
rz(1.596405) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.77869271) q[3];
sx q[3];
rz(-2.7540996) q[3];
sx q[3];
rz(-0.029822895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.80660194) q[2];
sx q[2];
rz(-1.0706173) q[2];
sx q[2];
rz(-0.09207329) q[2];
rz(-2.4798685) q[3];
sx q[3];
rz(-2.3501553) q[3];
sx q[3];
rz(1.7592643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6732366) q[0];
sx q[0];
rz(-1.3156923) q[0];
sx q[0];
rz(-3.1150505) q[0];
rz(-0.87310711) q[1];
sx q[1];
rz(-1.1353506) q[1];
sx q[1];
rz(2.81566) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8317141) q[0];
sx q[0];
rz(-1.1967812) q[0];
sx q[0];
rz(-0.563234) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8750149) q[2];
sx q[2];
rz(-1.1433257) q[2];
sx q[2];
rz(-1.3104591) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1205475) q[1];
sx q[1];
rz(-1.7056587) q[1];
sx q[1];
rz(0.30121505) q[1];
rz(-0.16485729) q[3];
sx q[3];
rz(-0.10247173) q[3];
sx q[3];
rz(-0.41302478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.59763336) q[2];
sx q[2];
rz(-1.8240857) q[2];
sx q[2];
rz(2.4874172) q[2];
rz(1.4298965) q[3];
sx q[3];
rz(-1.9734029) q[3];
sx q[3];
rz(-2.6005319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8687246) q[0];
sx q[0];
rz(-1.6690212) q[0];
sx q[0];
rz(2.4196999) q[0];
rz(-1.7294653) q[1];
sx q[1];
rz(-2.3528603) q[1];
sx q[1];
rz(3.022335) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28173791) q[0];
sx q[0];
rz(-0.87346948) q[0];
sx q[0];
rz(1.9182693) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.10995933) q[2];
sx q[2];
rz(-1.409515) q[2];
sx q[2];
rz(1.9732628) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.99265656) q[1];
sx q[1];
rz(-1.438237) q[1];
sx q[1];
rz(-2.0723144) q[1];
rz(2.4339606) q[3];
sx q[3];
rz(-1.283657) q[3];
sx q[3];
rz(-1.6361145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.44935903) q[2];
sx q[2];
rz(-1.8211775) q[2];
sx q[2];
rz(1.7822441) q[2];
rz(-2.3826777) q[3];
sx q[3];
rz(-0.24154285) q[3];
sx q[3];
rz(2.6045077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83157241) q[0];
sx q[0];
rz(-2.4825403) q[0];
sx q[0];
rz(-2.0781562) q[0];
rz(2.8670782) q[1];
sx q[1];
rz(-1.2083222) q[1];
sx q[1];
rz(2.2559821) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5078686) q[0];
sx q[0];
rz(-1.8025724) q[0];
sx q[0];
rz(-2.6404882) q[0];
rz(-pi) q[1];
rz(1.5468855) q[2];
sx q[2];
rz(-2.5619321) q[2];
sx q[2];
rz(-1.5469345) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1427666) q[1];
sx q[1];
rz(-0.66122675) q[1];
sx q[1];
rz(2.2426474) q[1];
rz(-pi) q[2];
rz(-2.4212491) q[3];
sx q[3];
rz(-0.90984905) q[3];
sx q[3];
rz(1.1870445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.72835913) q[2];
sx q[2];
rz(-0.76247549) q[2];
sx q[2];
rz(-1.1317066) q[2];
rz(-2.0570095) q[3];
sx q[3];
rz(-2.0621433) q[3];
sx q[3];
rz(-1.926698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30329147) q[0];
sx q[0];
rz(-1.6475995) q[0];
sx q[0];
rz(-1.0820748) q[0];
rz(1.8661631) q[1];
sx q[1];
rz(-2.137303) q[1];
sx q[1];
rz(-1.1358322) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71400434) q[0];
sx q[0];
rz(-1.7363318) q[0];
sx q[0];
rz(1.0323314) q[0];
rz(-pi) q[1];
rz(-0.047770569) q[2];
sx q[2];
rz(-2.7063745) q[2];
sx q[2];
rz(1.4620632) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.218704) q[1];
sx q[1];
rz(-1.6195546) q[1];
sx q[1];
rz(-2.7152275) q[1];
x q[2];
rz(-1.0335835) q[3];
sx q[3];
rz(-2.4164003) q[3];
sx q[3];
rz(2.4560526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.11848005) q[2];
sx q[2];
rz(-1.8818405) q[2];
sx q[2];
rz(0.87289587) q[2];
rz(-2.2980799) q[3];
sx q[3];
rz(-2.8819363) q[3];
sx q[3];
rz(-0.18994722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89235598) q[0];
sx q[0];
rz(-1.3441688) q[0];
sx q[0];
rz(1.2783485) q[0];
rz(1.0247914) q[1];
sx q[1];
rz(-2.0147851) q[1];
sx q[1];
rz(1.9445673) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9997864) q[0];
sx q[0];
rz(-2.1394661) q[0];
sx q[0];
rz(2.6151711) q[0];
x q[1];
rz(-3.0212101) q[2];
sx q[2];
rz(-1.8070081) q[2];
sx q[2];
rz(-1.2506968) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5843643) q[1];
sx q[1];
rz(-0.96712501) q[1];
sx q[1];
rz(-0.61324688) q[1];
rz(-pi) q[2];
x q[2];
rz(0.50935575) q[3];
sx q[3];
rz(-1.9037814) q[3];
sx q[3];
rz(-1.8978564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8355576) q[2];
sx q[2];
rz(-2.4287537) q[2];
sx q[2];
rz(-0.79997921) q[2];
rz(-1.9647313) q[3];
sx q[3];
rz(-1.7088944) q[3];
sx q[3];
rz(2.1879788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4326614) q[0];
sx q[0];
rz(-2.9928757) q[0];
sx q[0];
rz(-2.3401674) q[0];
rz(2.6196383) q[1];
sx q[1];
rz(-2.3028761) q[1];
sx q[1];
rz(0.16470673) q[1];
rz(1.9258826) q[2];
sx q[2];
rz(-1.2181031) q[2];
sx q[2];
rz(1.9545771) q[2];
rz(-2.342631) q[3];
sx q[3];
rz(-2.3286455) q[3];
sx q[3];
rz(0.13959985) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];