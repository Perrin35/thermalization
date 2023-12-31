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
rz(0.82984501) q[0];
rz(3.9217477) q[1];
sx q[1];
rz(5.2182066) q[1];
sx q[1];
rz(10.301104) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0759461) q[0];
sx q[0];
rz(-1.8008721) q[0];
sx q[0];
rz(-0.40212698) q[0];
x q[1];
rz(-2.1562188) q[2];
sx q[2];
rz(-1.4101763) q[2];
sx q[2];
rz(-2.0656245) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0383366) q[1];
sx q[1];
rz(-1.5557319) q[1];
sx q[1];
rz(-1.2407833) q[1];
x q[2];
rz(2.797804) q[3];
sx q[3];
rz(-0.4292092) q[3];
sx q[3];
rz(2.6478298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9709388) q[2];
sx q[2];
rz(-1.8654414) q[2];
sx q[2];
rz(-0.7286287) q[2];
rz(2.6206) q[3];
sx q[3];
rz(-0.96121585) q[3];
sx q[3];
rz(0.20761028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8347297) q[0];
sx q[0];
rz(-1.1704209) q[0];
sx q[0];
rz(1.1215425) q[0];
rz(2.8858378) q[1];
sx q[1];
rz(-1.6633727) q[1];
sx q[1];
rz(2.2671525) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53423184) q[0];
sx q[0];
rz(-0.53593862) q[0];
sx q[0];
rz(2.1952654) q[0];
rz(0.15408709) q[2];
sx q[2];
rz(-1.0006957) q[2];
sx q[2];
rz(-2.4947583) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.16972152) q[1];
sx q[1];
rz(-2.019878) q[1];
sx q[1];
rz(-0.68193087) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8401243) q[3];
sx q[3];
rz(-1.519031) q[3];
sx q[3];
rz(-0.89382899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4008537) q[2];
sx q[2];
rz(-0.48745552) q[2];
sx q[2];
rz(-2.7056616) q[2];
rz(0.68108264) q[3];
sx q[3];
rz(-2.3705132) q[3];
sx q[3];
rz(-2.7387103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(-0.23713672) q[0];
sx q[0];
rz(-0.28770068) q[0];
sx q[0];
rz(1.0748192) q[0];
rz(-2.3020321) q[1];
sx q[1];
rz(-0.81937516) q[1];
sx q[1];
rz(-2.7456465) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9112644) q[0];
sx q[0];
rz(-1.0697782) q[0];
sx q[0];
rz(2.688174) q[0];
x q[1];
rz(1.5244353) q[2];
sx q[2];
rz(-2.8472387) q[2];
sx q[2];
rz(-0.72325318) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1173276) q[1];
sx q[1];
rz(-1.5803442) q[1];
sx q[1];
rz(-2.4584103) q[1];
rz(-pi) q[2];
rz(2.8860502) q[3];
sx q[3];
rz(-1.5934172) q[3];
sx q[3];
rz(2.8850151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.6039156) q[2];
sx q[2];
rz(-1.2664412) q[2];
sx q[2];
rz(0.23920693) q[2];
rz(3.0662597) q[3];
sx q[3];
rz(-2.0276666) q[3];
sx q[3];
rz(1.8384365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72162119) q[0];
sx q[0];
rz(-0.8525089) q[0];
sx q[0];
rz(2.3216632) q[0];
rz(-2.6539102) q[1];
sx q[1];
rz(-0.90355021) q[1];
sx q[1];
rz(2.908169) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48196402) q[0];
sx q[0];
rz(-3.0273962) q[0];
sx q[0];
rz(-1.0114848) q[0];
x q[1];
rz(0.64220631) q[2];
sx q[2];
rz(-1.7608479) q[2];
sx q[2];
rz(1.9501291) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.1366795) q[1];
sx q[1];
rz(-2.3466952) q[1];
sx q[1];
rz(-0.48308785) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.78419533) q[3];
sx q[3];
rz(-1.2894221) q[3];
sx q[3];
rz(-0.97660645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0908115) q[2];
sx q[2];
rz(-1.2819141) q[2];
sx q[2];
rz(-0.74742571) q[2];
rz(-2.9181972) q[3];
sx q[3];
rz(-2.5441393) q[3];
sx q[3];
rz(-2.8994765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4500047) q[0];
sx q[0];
rz(-1.4331899) q[0];
sx q[0];
rz(-3.0773556) q[0];
rz(-0.94379395) q[1];
sx q[1];
rz(-0.72729021) q[1];
sx q[1];
rz(0.76104004) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8781389) q[0];
sx q[0];
rz(-1.8315151) q[0];
sx q[0];
rz(0.031147416) q[0];
x q[1];
rz(2.9391187) q[2];
sx q[2];
rz(-1.7624117) q[2];
sx q[2];
rz(1.7970049) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0423454) q[1];
sx q[1];
rz(-1.595101) q[1];
sx q[1];
rz(-2.821032) q[1];
rz(-pi) q[2];
rz(-2.8588572) q[3];
sx q[3];
rz(-1.3021819) q[3];
sx q[3];
rz(-2.3408567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3349907) q[2];
sx q[2];
rz(-2.0709753) q[2];
sx q[2];
rz(-3.0495194) q[2];
rz(-2.4798685) q[3];
sx q[3];
rz(-0.79143733) q[3];
sx q[3];
rz(1.3823284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6732366) q[0];
sx q[0];
rz(-1.8259003) q[0];
sx q[0];
rz(-0.026542149) q[0];
rz(-0.87310711) q[1];
sx q[1];
rz(-2.006242) q[1];
sx q[1];
rz(0.32593265) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8317141) q[0];
sx q[0];
rz(-1.9448115) q[0];
sx q[0];
rz(0.563234) q[0];
rz(-pi) q[1];
x q[1];
rz(0.44547703) q[2];
sx q[2];
rz(-1.8468841) q[2];
sx q[2];
rz(2.7518227) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.50800486) q[1];
sx q[1];
rz(-1.8691917) q[1];
sx q[1];
rz(-1.7119346) q[1];
rz(-pi) q[2];
rz(1.5876706) q[3];
sx q[3];
rz(-1.4697187) q[3];
sx q[3];
rz(-0.24731393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.59763336) q[2];
sx q[2];
rz(-1.3175069) q[2];
sx q[2];
rz(-0.65417543) q[2];
rz(-1.7116961) q[3];
sx q[3];
rz(-1.9734029) q[3];
sx q[3];
rz(0.54106075) q[3];
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
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27286801) q[0];
sx q[0];
rz(-1.4725715) q[0];
sx q[0];
rz(-0.72189271) q[0];
rz(1.4121274) q[1];
sx q[1];
rz(-0.78873235) q[1];
sx q[1];
rz(-3.022335) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0810453) q[0];
sx q[0];
rz(-1.8348872) q[0];
sx q[0];
rz(2.413785) q[0];
rz(-1.4085521) q[2];
sx q[2];
rz(-1.4622697) q[2];
sx q[2];
rz(0.42019368) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6357928) q[1];
sx q[1];
rz(-2.0675106) q[1];
sx q[1];
rz(-0.15091166) q[1];
rz(-pi) q[2];
rz(1.2001541) q[3];
sx q[3];
rz(-0.89768411) q[3];
sx q[3];
rz(2.9692269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6922336) q[2];
sx q[2];
rz(-1.3204152) q[2];
sx q[2];
rz(1.7822441) q[2];
rz(-0.75891495) q[3];
sx q[3];
rz(-2.9000498) q[3];
sx q[3];
rz(2.6045077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.3100202) q[0];
sx q[0];
rz(-2.4825403) q[0];
sx q[0];
rz(-1.0634364) q[0];
rz(0.27451441) q[1];
sx q[1];
rz(-1.9332705) q[1];
sx q[1];
rz(2.2559821) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63372406) q[0];
sx q[0];
rz(-1.8025724) q[0];
sx q[0];
rz(-0.50110441) q[0];
x q[1];
rz(0.015651264) q[2];
sx q[2];
rz(-0.99132292) q[2];
sx q[2];
rz(1.5755115) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9319716) q[1];
sx q[1];
rz(-2.0721657) q[1];
sx q[1];
rz(0.45100905) q[1];
rz(-pi) q[2];
x q[2];
rz(0.7203436) q[3];
sx q[3];
rz(-2.2317436) q[3];
sx q[3];
rz(-1.1870445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.72835913) q[2];
sx q[2];
rz(-0.76247549) q[2];
sx q[2];
rz(1.1317066) q[2];
rz(2.0570095) q[3];
sx q[3];
rz(-2.0621433) q[3];
sx q[3];
rz(-1.2148946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8383012) q[0];
sx q[0];
rz(-1.6475995) q[0];
sx q[0];
rz(-2.0595179) q[0];
rz(-1.2754296) q[1];
sx q[1];
rz(-2.137303) q[1];
sx q[1];
rz(-1.1358322) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4275883) q[0];
sx q[0];
rz(-1.7363318) q[0];
sx q[0];
rz(2.1092613) q[0];
rz(-2.7068107) q[2];
sx q[2];
rz(-1.5909305) q[2];
sx q[2];
rz(0.15205631) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3867934) q[1];
sx q[1];
rz(-2.7126185) q[1];
sx q[1];
rz(3.0241443) q[1];
rz(2.2215861) q[3];
sx q[3];
rz(-1.2244867) q[3];
sx q[3];
rz(-1.3045834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0231126) q[2];
sx q[2];
rz(-1.8818405) q[2];
sx q[2];
rz(-0.87289587) q[2];
rz(-0.84351271) q[3];
sx q[3];
rz(-0.25965634) q[3];
sx q[3];
rz(-0.18994722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2492367) q[0];
sx q[0];
rz(-1.3441688) q[0];
sx q[0];
rz(-1.8632442) q[0];
rz(1.0247914) q[1];
sx q[1];
rz(-1.1268076) q[1];
sx q[1];
rz(-1.9445673) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68162936) q[0];
sx q[0];
rz(-2.3869793) q[0];
sx q[0];
rz(0.90453903) q[0];
rz(-1.808666) q[2];
sx q[2];
rz(-1.6878205) q[2];
sx q[2];
rz(0.34840096) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.69233209) q[1];
sx q[1];
rz(-0.83220607) q[1];
sx q[1];
rz(2.2663121) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6322369) q[3];
sx q[3];
rz(-1.2378113) q[3];
sx q[3];
rz(-1.8978564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8355576) q[2];
sx q[2];
rz(-0.71283895) q[2];
sx q[2];
rz(-2.3416134) q[2];
rz(1.1768613) q[3];
sx q[3];
rz(-1.4326982) q[3];
sx q[3];
rz(0.95361382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4326614) q[0];
sx q[0];
rz(-0.14871696) q[0];
sx q[0];
rz(0.8014252) q[0];
rz(0.52195436) q[1];
sx q[1];
rz(-0.83871651) q[1];
sx q[1];
rz(-2.9768859) q[1];
rz(2.7675046) q[2];
sx q[2];
rz(-1.2384407) q[2];
sx q[2];
rz(0.25638914) q[2];
rz(-2.5064777) q[3];
sx q[3];
rz(-2.1182346) q[3];
sx q[3];
rz(2.3253141) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
