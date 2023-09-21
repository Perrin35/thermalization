OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.61385566) q[0];
sx q[0];
rz(-1.6439438) q[0];
sx q[0];
rz(-0.82984501) q[0];
rz(0.78015503) q[1];
sx q[1];
rz(-2.0766139) q[1];
sx q[1];
rz(-0.87632626) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5397545) q[0];
sx q[0];
rz(-1.9617426) q[0];
sx q[0];
rz(-1.3215617) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8560156) q[2];
sx q[2];
rz(-2.5370295) q[2];
sx q[2];
rz(2.4100458) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6688924) q[1];
sx q[1];
rz(-1.9007705) q[1];
sx q[1];
rz(-0.015923576) q[1];
rz(-pi) q[2];
rz(-2.734749) q[3];
sx q[3];
rz(-1.4300656) q[3];
sx q[3];
rz(-1.391747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.17065389) q[2];
sx q[2];
rz(-1.2761513) q[2];
sx q[2];
rz(0.7286287) q[2];
rz(-2.6206) q[3];
sx q[3];
rz(-2.1803768) q[3];
sx q[3];
rz(-2.9339824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.306863) q[0];
sx q[0];
rz(-1.1704209) q[0];
sx q[0];
rz(-1.1215425) q[0];
rz(-2.8858378) q[1];
sx q[1];
rz(-1.6633727) q[1];
sx q[1];
rz(0.87444011) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.550299) q[0];
sx q[0];
rz(-1.2676139) q[0];
sx q[0];
rz(-1.1217872) q[0];
rz(-pi) q[1];
x q[1];
rz(0.15408709) q[2];
sx q[2];
rz(-2.140897) q[2];
sx q[2];
rz(-0.64683435) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9718711) q[1];
sx q[1];
rz(-2.019878) q[1];
sx q[1];
rz(0.68193087) q[1];
x q[2];
rz(3.087895) q[3];
sx q[3];
rz(-1.8397545) q[3];
sx q[3];
rz(0.69124903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.740739) q[2];
sx q[2];
rz(-0.48745552) q[2];
sx q[2];
rz(2.7056616) q[2];
rz(0.68108264) q[3];
sx q[3];
rz(-0.77107945) q[3];
sx q[3];
rz(-0.40288231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23713672) q[0];
sx q[0];
rz(-0.28770068) q[0];
sx q[0];
rz(2.0667734) q[0];
rz(-2.3020321) q[1];
sx q[1];
rz(-2.3222175) q[1];
sx q[1];
rz(-0.39594617) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0229189) q[0];
sx q[0];
rz(-0.6624822) q[0];
sx q[0];
rz(-0.89612095) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.01404889) q[2];
sx q[2];
rz(-1.2767681) q[2];
sx q[2];
rz(-0.77169466) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1173276) q[1];
sx q[1];
rz(-1.5612484) q[1];
sx q[1];
rz(0.68318232) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8860502) q[3];
sx q[3];
rz(-1.5481755) q[3];
sx q[3];
rz(-0.25657755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6039156) q[2];
sx q[2];
rz(-1.8751514) q[2];
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
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4199715) q[0];
sx q[0];
rz(-2.2890838) q[0];
sx q[0];
rz(-2.3216632) q[0];
rz(2.6539102) q[1];
sx q[1];
rz(-0.90355021) q[1];
sx q[1];
rz(-2.908169) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.061302) q[0];
sx q[0];
rz(-1.6675321) q[0];
sx q[0];
rz(3.0808099) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8065622) q[2];
sx q[2];
rz(-2.1996017) q[2];
sx q[2];
rz(2.6218888) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.1366795) q[1];
sx q[1];
rz(-2.3466952) q[1];
sx q[1];
rz(-0.48308785) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1831746) q[3];
sx q[3];
rz(-0.82510199) q[3];
sx q[3];
rz(-0.86442282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0507811) q[2];
sx q[2];
rz(-1.2819141) q[2];
sx q[2];
rz(2.3941669) q[2];
rz(2.9181972) q[3];
sx q[3];
rz(-0.59745336) q[3];
sx q[3];
rz(-2.8994765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(-1.4500047) q[0];
sx q[0];
rz(-1.7084028) q[0];
sx q[0];
rz(-0.064237021) q[0];
rz(2.1977987) q[1];
sx q[1];
rz(-2.4143024) q[1];
sx q[1];
rz(-0.76104004) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1431664) q[0];
sx q[0];
rz(-0.26253065) q[0];
sx q[0];
rz(1.6869998) q[0];
x q[1];
rz(-1.3752851) q[2];
sx q[2];
rz(-1.7695145) q[2];
sx q[2];
rz(-2.9544601) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4796175) q[1];
sx q[1];
rz(-1.2503337) q[1];
sx q[1];
rz(1.5451876) q[1];
x q[2];
rz(-1.8499591) q[3];
sx q[3];
rz(-1.2984635) q[3];
sx q[3];
rz(-0.84701049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.80660194) q[2];
sx q[2];
rz(-1.0706173) q[2];
sx q[2];
rz(0.09207329) q[2];
rz(0.66172415) q[3];
sx q[3];
rz(-2.3501553) q[3];
sx q[3];
rz(1.7592643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46835607) q[0];
sx q[0];
rz(-1.3156923) q[0];
sx q[0];
rz(-0.026542149) q[0];
rz(-2.2684855) q[1];
sx q[1];
rz(-2.006242) q[1];
sx q[1];
rz(-0.32593265) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3098785) q[0];
sx q[0];
rz(-1.1967812) q[0];
sx q[0];
rz(0.563234) q[0];
rz(-pi) q[1];
rz(-0.58165254) q[2];
sx q[2];
rz(-0.51917167) q[2];
sx q[2];
rz(0.66228629) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0210452) q[1];
sx q[1];
rz(-1.435934) q[1];
sx q[1];
rz(2.8403776) q[1];
x q[2];
rz(-0.10109191) q[3];
sx q[3];
rz(-1.5875845) q[3];
sx q[3];
rz(1.8198131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5439593) q[2];
sx q[2];
rz(-1.8240857) q[2];
sx q[2];
rz(0.65417543) q[2];
rz(-1.4298965) q[3];
sx q[3];
rz(-1.9734029) q[3];
sx q[3];
rz(2.6005319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8687246) q[0];
sx q[0];
rz(-1.4725715) q[0];
sx q[0];
rz(0.72189271) q[0];
rz(-1.7294653) q[1];
sx q[1];
rz(-2.3528603) q[1];
sx q[1];
rz(3.022335) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0605474) q[0];
sx q[0];
rz(-1.8348872) q[0];
sx q[0];
rz(-0.72780769) q[0];
rz(-1.7330405) q[2];
sx q[2];
rz(-1.4622697) q[2];
sx q[2];
rz(2.721399) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3268986) q[1];
sx q[1];
rz(-2.6242995) q[1];
sx q[1];
rz(1.8413522) q[1];
rz(1.9414385) q[3];
sx q[3];
rz(-0.89768411) q[3];
sx q[3];
rz(0.17236575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6922336) q[2];
sx q[2];
rz(-1.3204152) q[2];
sx q[2];
rz(-1.3593486) q[2];
rz(-2.3826777) q[3];
sx q[3];
rz(-2.9000498) q[3];
sx q[3];
rz(0.53708491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3100202) q[0];
sx q[0];
rz(-0.65905237) q[0];
sx q[0];
rz(1.0634364) q[0];
rz(-0.27451441) q[1];
sx q[1];
rz(-1.2083222) q[1];
sx q[1];
rz(2.2559821) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5078686) q[0];
sx q[0];
rz(-1.8025724) q[0];
sx q[0];
rz(0.50110441) q[0];
rz(1.5468855) q[2];
sx q[2];
rz(-2.5619321) q[2];
sx q[2];
rz(1.5946582) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0091128) q[1];
sx q[1];
rz(-1.9630034) q[1];
sx q[1];
rz(1.023804) q[1];
x q[2];
rz(-0.86730154) q[3];
sx q[3];
rz(-2.2059545) q[3];
sx q[3];
rz(0.22658539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.72835913) q[2];
sx q[2];
rz(-0.76247549) q[2];
sx q[2];
rz(1.1317066) q[2];
rz(-2.0570095) q[3];
sx q[3];
rz(-1.0794493) q[3];
sx q[3];
rz(-1.2148946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30329147) q[0];
sx q[0];
rz(-1.4939932) q[0];
sx q[0];
rz(1.0820748) q[0];
rz(-1.8661631) q[1];
sx q[1];
rz(-1.0042896) q[1];
sx q[1];
rz(-1.1358322) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58763114) q[0];
sx q[0];
rz(-2.580664) q[0];
sx q[0];
rz(1.8857303) q[0];
rz(-2.7068107) q[2];
sx q[2];
rz(-1.5506622) q[2];
sx q[2];
rz(2.9895363) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.62577265) q[1];
sx q[1];
rz(-1.9966218) q[1];
sx q[1];
rz(-1.6243402) q[1];
x q[2];
rz(-2.7157852) q[3];
sx q[3];
rz(-2.1771181) q[3];
sx q[3];
rz(-3.1283034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0231126) q[2];
sx q[2];
rz(-1.8818405) q[2];
sx q[2];
rz(2.2686968) q[2];
rz(2.2980799) q[3];
sx q[3];
rz(-0.25965634) q[3];
sx q[3];
rz(2.9516454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89235598) q[0];
sx q[0];
rz(-1.7974239) q[0];
sx q[0];
rz(-1.2783485) q[0];
rz(-1.0247914) q[1];
sx q[1];
rz(-1.1268076) q[1];
sx q[1];
rz(1.9445673) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14180627) q[0];
sx q[0];
rz(-2.1394661) q[0];
sx q[0];
rz(0.52642157) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1080152) q[2];
sx q[2];
rz(-0.26460755) q[2];
sx q[2];
rz(-0.77361425) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7750818) q[1];
sx q[1];
rz(-2.0644036) q[1];
sx q[1];
rz(0.8702741) q[1];
rz(-pi) q[2];
rz(-1.9479806) q[3];
sx q[3];
rz(-1.0918655) q[3];
sx q[3];
rz(0.14648986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8355576) q[2];
sx q[2];
rz(-2.4287537) q[2];
sx q[2];
rz(-0.79997921) q[2];
rz(1.9647313) q[3];
sx q[3];
rz(-1.4326982) q[3];
sx q[3];
rz(-0.95361382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70893127) q[0];
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
rz(2.5064777) q[3];
sx q[3];
rz(-1.0233581) q[3];
sx q[3];
rz(-0.81627853) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
