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
rz(0.78015503) q[1];
sx q[1];
rz(-2.0766139) q[1];
sx q[1];
rz(2.2652664) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5397545) q[0];
sx q[0];
rz(-1.1798501) q[0];
sx q[0];
rz(-1.3215617) q[0];
rz(-pi) q[1];
rz(-1.285577) q[2];
sx q[2];
rz(-2.5370295) q[2];
sx q[2];
rz(2.4100458) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7180011) q[1];
sx q[1];
rz(-0.33034409) q[1];
sx q[1];
rz(1.6172536) q[1];
x q[2];
rz(1.7238486) q[3];
sx q[3];
rz(-1.1682086) q[3];
sx q[3];
rz(3.0229085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9709388) q[2];
sx q[2];
rz(-1.8654414) q[2];
sx q[2];
rz(2.412964) q[2];
rz(2.6206) q[3];
sx q[3];
rz(-2.1803768) q[3];
sx q[3];
rz(-0.20761028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8347297) q[0];
sx q[0];
rz(-1.1704209) q[0];
sx q[0];
rz(-1.1215425) q[0];
rz(-0.25575486) q[1];
sx q[1];
rz(-1.47822) q[1];
sx q[1];
rz(0.87444011) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.550299) q[0];
sx q[0];
rz(-1.2676139) q[0];
sx q[0];
rz(2.0198054) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3358243) q[2];
sx q[2];
rz(-0.58832303) q[2];
sx q[2];
rz(-2.7749643) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2494617) q[1];
sx q[1];
rz(-2.3453237) q[1];
sx q[1];
rz(-0.65278058) q[1];
rz(0.053697649) q[3];
sx q[3];
rz(-1.3018381) q[3];
sx q[3];
rz(-2.4503436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.740739) q[2];
sx q[2];
rz(-2.6541371) q[2];
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
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23713672) q[0];
sx q[0];
rz(-2.853892) q[0];
sx q[0];
rz(2.0667734) q[0];
rz(-0.83956051) q[1];
sx q[1];
rz(-0.81937516) q[1];
sx q[1];
rz(-0.39594617) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9112644) q[0];
sx q[0];
rz(-2.0718144) q[0];
sx q[0];
rz(0.45341861) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1275438) q[2];
sx q[2];
rz(-1.2767681) q[2];
sx q[2];
rz(-2.369898) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.695895) q[1];
sx q[1];
rz(-2.2539415) q[1];
sx q[1];
rz(1.5584857) q[1];
x q[2];
rz(0.08926908) q[3];
sx q[3];
rz(-0.25651989) q[3];
sx q[3];
rz(-1.741011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5376771) q[2];
sx q[2];
rz(-1.2664412) q[2];
sx q[2];
rz(-2.9023857) q[2];
rz(3.0662597) q[3];
sx q[3];
rz(-2.0276666) q[3];
sx q[3];
rz(1.8384365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4199715) q[0];
sx q[0];
rz(-0.8525089) q[0];
sx q[0];
rz(-0.81992942) q[0];
rz(0.48768249) q[1];
sx q[1];
rz(-2.2380424) q[1];
sx q[1];
rz(-2.908169) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48196402) q[0];
sx q[0];
rz(-3.0273962) q[0];
sx q[0];
rz(1.0114848) q[0];
rz(-pi) q[1];
rz(-0.31077023) q[2];
sx q[2];
rz(-0.66590532) q[2];
sx q[2];
rz(0.13194612) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.0049131752) q[1];
sx q[1];
rz(-0.79489743) q[1];
sx q[1];
rz(-2.6585048) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7531284) q[3];
sx q[3];
rz(-0.82287517) q[3];
sx q[3];
rz(-0.32271656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0908115) q[2];
sx q[2];
rz(-1.2819141) q[2];
sx q[2];
rz(0.74742571) q[2];
rz(-0.22339544) q[3];
sx q[3];
rz(-0.59745336) q[3];
sx q[3];
rz(0.24211611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4500047) q[0];
sx q[0];
rz(-1.7084028) q[0];
sx q[0];
rz(3.0773556) q[0];
rz(2.1977987) q[1];
sx q[1];
rz(-2.4143024) q[1];
sx q[1];
rz(2.3805526) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9984263) q[0];
sx q[0];
rz(-0.26253065) q[0];
sx q[0];
rz(-1.4545928) q[0];
x q[1];
rz(0.76743482) q[2];
sx q[2];
rz(-0.27786294) q[2];
sx q[2];
rz(-2.1674736) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.099247301) q[1];
sx q[1];
rz(-1.5464916) q[1];
sx q[1];
rz(-2.821032) q[1];
rz(-pi) q[2];
rz(-2.8588572) q[3];
sx q[3];
rz(-1.8394107) q[3];
sx q[3];
rz(2.3408567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3349907) q[2];
sx q[2];
rz(-1.0706173) q[2];
sx q[2];
rz(-3.0495194) q[2];
rz(-0.66172415) q[3];
sx q[3];
rz(-0.79143733) q[3];
sx q[3];
rz(-1.3823284) q[3];
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
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46835607) q[0];
sx q[0];
rz(-1.3156923) q[0];
sx q[0];
rz(0.026542149) q[0];
rz(2.2684855) q[1];
sx q[1];
rz(-1.1353506) q[1];
sx q[1];
rz(2.81566) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3098785) q[0];
sx q[0];
rz(-1.1967812) q[0];
sx q[0];
rz(2.5783587) q[0];
rz(-2.6961156) q[2];
sx q[2];
rz(-1.8468841) q[2];
sx q[2];
rz(-0.38976994) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.50800486) q[1];
sx q[1];
rz(-1.8691917) q[1];
sx q[1];
rz(-1.7119346) q[1];
rz(-pi) q[2];
x q[2];
rz(0.16485729) q[3];
sx q[3];
rz(-0.10247173) q[3];
sx q[3];
rz(0.41302478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.59763336) q[2];
sx q[2];
rz(-1.8240857) q[2];
sx q[2];
rz(-2.4874172) q[2];
rz(1.7116961) q[3];
sx q[3];
rz(-1.9734029) q[3];
sx q[3];
rz(-0.54106075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27286801) q[0];
sx q[0];
rz(-1.4725715) q[0];
sx q[0];
rz(-2.4196999) q[0];
rz(-1.4121274) q[1];
sx q[1];
rz(-2.3528603) q[1];
sx q[1];
rz(0.11925764) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3463466) q[0];
sx q[0];
rz(-0.7659142) q[0];
sx q[0];
rz(-0.38608293) q[0];
rz(-2.1642045) q[2];
sx q[2];
rz(-0.19492976) q[2];
sx q[2];
rz(0.56602636) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6357928) q[1];
sx q[1];
rz(-1.074082) q[1];
sx q[1];
rz(2.990681) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.70763208) q[3];
sx q[3];
rz(-1.283657) q[3];
sx q[3];
rz(-1.6361145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.44935903) q[2];
sx q[2];
rz(-1.3204152) q[2];
sx q[2];
rz(-1.7822441) q[2];
rz(0.75891495) q[3];
sx q[3];
rz(-0.24154285) q[3];
sx q[3];
rz(-0.53708491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3100202) q[0];
sx q[0];
rz(-0.65905237) q[0];
sx q[0];
rz(-1.0634364) q[0];
rz(-0.27451441) q[1];
sx q[1];
rz(-1.2083222) q[1];
sx q[1];
rz(-0.88561052) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5078686) q[0];
sx q[0];
rz(-1.8025724) q[0];
sx q[0];
rz(-0.50110441) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5947072) q[2];
sx q[2];
rz(-0.57966053) q[2];
sx q[2];
rz(1.5469345) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.209621) q[1];
sx q[1];
rz(-1.069427) q[1];
sx q[1];
rz(-0.45100905) q[1];
x q[2];
rz(-0.76836821) q[3];
sx q[3];
rz(-1.0230912) q[3];
sx q[3];
rz(0.87793575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.72835913) q[2];
sx q[2];
rz(-0.76247549) q[2];
sx q[2];
rz(1.1317066) q[2];
rz(-1.0845832) q[3];
sx q[3];
rz(-2.0621433) q[3];
sx q[3];
rz(1.926698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30329147) q[0];
sx q[0];
rz(-1.4939932) q[0];
sx q[0];
rz(-1.0820748) q[0];
rz(-1.2754296) q[1];
sx q[1];
rz(-2.137303) q[1];
sx q[1];
rz(2.0057604) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5539615) q[0];
sx q[0];
rz(-0.56092867) q[0];
sx q[0];
rz(1.8857303) q[0];
x q[1];
rz(-3.0938221) q[2];
sx q[2];
rz(-0.43521817) q[2];
sx q[2];
rz(-1.6795295) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.75479924) q[1];
sx q[1];
rz(-0.42897412) q[1];
sx q[1];
rz(-3.0241443) q[1];
x q[2];
rz(-2.2215861) q[3];
sx q[3];
rz(-1.917106) q[3];
sx q[3];
rz(1.8370093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0231126) q[2];
sx q[2];
rz(-1.2597522) q[2];
sx q[2];
rz(0.87289587) q[2];
rz(-0.84351271) q[3];
sx q[3];
rz(-2.8819363) q[3];
sx q[3];
rz(-2.9516454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89235598) q[0];
sx q[0];
rz(-1.3441688) q[0];
sx q[0];
rz(-1.8632442) q[0];
rz(2.1168013) q[1];
sx q[1];
rz(-1.1268076) q[1];
sx q[1];
rz(1.9445673) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4599633) q[0];
sx q[0];
rz(-2.3869793) q[0];
sx q[0];
rz(2.2370536) q[0];
x q[1];
rz(1.3329266) q[2];
sx q[2];
rz(-1.6878205) q[2];
sx q[2];
rz(-2.7931917) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5572284) q[1];
sx q[1];
rz(-2.1744676) q[1];
sx q[1];
rz(0.61324688) q[1];
rz(-pi) q[2];
x q[2];
rz(0.61693807) q[3];
sx q[3];
rz(-2.5411798) q[3];
sx q[3];
rz(2.2850349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3060351) q[2];
sx q[2];
rz(-2.4287537) q[2];
sx q[2];
rz(2.3416134) q[2];
rz(1.1768613) q[3];
sx q[3];
rz(-1.7088944) q[3];
sx q[3];
rz(2.1879788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4326614) q[0];
sx q[0];
rz(-2.9928757) q[0];
sx q[0];
rz(-2.3401674) q[0];
rz(-2.6196383) q[1];
sx q[1];
rz(-0.83871651) q[1];
sx q[1];
rz(-2.9768859) q[1];
rz(2.3847053) q[2];
sx q[2];
rz(-0.49513985) q[2];
sx q[2];
rz(-2.0078299) q[2];
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