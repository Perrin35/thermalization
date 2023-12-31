OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.7735908) q[0];
sx q[0];
rz(-2.350783) q[0];
sx q[0];
rz(2.8074582) q[0];
rz(2.6842527) q[1];
sx q[1];
rz(-2.1973124) q[1];
sx q[1];
rz(-1.2184719) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2479808) q[0];
sx q[0];
rz(-0.75151822) q[0];
sx q[0];
rz(3.0889838) q[0];
rz(-pi) q[1];
rz(1.5300418) q[2];
sx q[2];
rz(-2.0353122) q[2];
sx q[2];
rz(-2.0704616) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6005046) q[1];
sx q[1];
rz(-1.0827912) q[1];
sx q[1];
rz(1.0469251) q[1];
x q[2];
rz(0.18913194) q[3];
sx q[3];
rz(-1.7889708) q[3];
sx q[3];
rz(-1.3355586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5228287) q[2];
sx q[2];
rz(-0.4814119) q[2];
sx q[2];
rz(0.5775601) q[2];
rz(1.1497568) q[3];
sx q[3];
rz(-1.7532319) q[3];
sx q[3];
rz(2.4770588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98786551) q[0];
sx q[0];
rz(-0.58652121) q[0];
sx q[0];
rz(-0.38744774) q[0];
rz(2.2024343) q[1];
sx q[1];
rz(-0.99717957) q[1];
sx q[1];
rz(1.739025) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99968602) q[0];
sx q[0];
rz(-1.5667331) q[0];
sx q[0];
rz(0.02948972) q[0];
rz(-pi) q[1];
rz(2.8005373) q[2];
sx q[2];
rz(-2.555072) q[2];
sx q[2];
rz(-1.7797433) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.1075588) q[1];
sx q[1];
rz(-0.65980136) q[1];
sx q[1];
rz(-2.4504513) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0121147) q[3];
sx q[3];
rz(-0.63614142) q[3];
sx q[3];
rz(-0.57487088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7188321) q[2];
sx q[2];
rz(-1.7451124) q[2];
sx q[2];
rz(2.823901) q[2];
rz(0.20673949) q[3];
sx q[3];
rz(-0.59967774) q[3];
sx q[3];
rz(0.81682214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-1.3690255) q[0];
sx q[0];
rz(-1.439753) q[0];
sx q[0];
rz(-1.7279708) q[0];
rz(-2.6638022) q[1];
sx q[1];
rz(-1.3505892) q[1];
sx q[1];
rz(-2.7405222) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7558407) q[0];
sx q[0];
rz(-2.7999561) q[0];
sx q[0];
rz(-1.2582448) q[0];
x q[1];
rz(0.6443278) q[2];
sx q[2];
rz(-1.4187078) q[2];
sx q[2];
rz(2.6459141) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9123147) q[1];
sx q[1];
rz(-0.73017987) q[1];
sx q[1];
rz(0.0095403949) q[1];
x q[2];
rz(1.1242261) q[3];
sx q[3];
rz(-0.43101573) q[3];
sx q[3];
rz(2.9197846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.59427375) q[2];
sx q[2];
rz(-1.5218364) q[2];
sx q[2];
rz(0.55580124) q[2];
rz(-2.1650971) q[3];
sx q[3];
rz(-0.54978168) q[3];
sx q[3];
rz(-0.78021375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7926086) q[0];
sx q[0];
rz(-1.5190834) q[0];
sx q[0];
rz(1.4439616) q[0];
rz(1.6216888) q[1];
sx q[1];
rz(-0.65602055) q[1];
sx q[1];
rz(2.8881883) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4454173) q[0];
sx q[0];
rz(-1.6735958) q[0];
sx q[0];
rz(-2.1131383) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1242712) q[2];
sx q[2];
rz(-1.0850731) q[2];
sx q[2];
rz(0.7893562) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4761915) q[1];
sx q[1];
rz(-0.63648495) q[1];
sx q[1];
rz(-0.16237662) q[1];
rz(2.4921791) q[3];
sx q[3];
rz(-1.489349) q[3];
sx q[3];
rz(-2.7922975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.110934) q[2];
sx q[2];
rz(-1.3867644) q[2];
sx q[2];
rz(-0.78732642) q[2];
rz(-2.2287255) q[3];
sx q[3];
rz(-0.74936167) q[3];
sx q[3];
rz(-2.1319938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.426429) q[0];
sx q[0];
rz(-0.63868317) q[0];
sx q[0];
rz(-0.062967904) q[0];
rz(-3.0175623) q[1];
sx q[1];
rz(-0.80563671) q[1];
sx q[1];
rz(2.6834992) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8040745) q[0];
sx q[0];
rz(-2.6989557) q[0];
sx q[0];
rz(-0.0054365693) q[0];
rz(-pi) q[1];
rz(-2.1796218) q[2];
sx q[2];
rz(-1.95032) q[2];
sx q[2];
rz(-0.48987197) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3921567) q[1];
sx q[1];
rz(-1.4538308) q[1];
sx q[1];
rz(2.3163296) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0714508) q[3];
sx q[3];
rz(-1.6541012) q[3];
sx q[3];
rz(1.7617776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9690341) q[2];
sx q[2];
rz(-0.92323747) q[2];
sx q[2];
rz(-0.22949533) q[2];
rz(3.138792) q[3];
sx q[3];
rz(-2.2715748) q[3];
sx q[3];
rz(-1.3389448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5795508) q[0];
sx q[0];
rz(-0.27844772) q[0];
sx q[0];
rz(2.2221185) q[0];
rz(3.0793076) q[1];
sx q[1];
rz(-2.1376164) q[1];
sx q[1];
rz(1.8744291) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29612449) q[0];
sx q[0];
rz(-2.2757029) q[0];
sx q[0];
rz(-1.0494997) q[0];
x q[1];
rz(-0.74395545) q[2];
sx q[2];
rz(-1.6774872) q[2];
sx q[2];
rz(-0.26847408) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8333203) q[1];
sx q[1];
rz(-2.7566524) q[1];
sx q[1];
rz(-1.0233378) q[1];
rz(-2.0503067) q[3];
sx q[3];
rz(-1.3538085) q[3];
sx q[3];
rz(-0.12245164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1798114) q[2];
sx q[2];
rz(-2.2634025) q[2];
sx q[2];
rz(0.58376694) q[2];
rz(0.70872712) q[3];
sx q[3];
rz(-1.3112336) q[3];
sx q[3];
rz(0.023199737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0631183) q[0];
sx q[0];
rz(-2.6874976) q[0];
sx q[0];
rz(-1.0725347) q[0];
rz(2.5947) q[1];
sx q[1];
rz(-1.2439589) q[1];
sx q[1];
rz(2.0297208) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3277153) q[0];
sx q[0];
rz(-1.0881256) q[0];
sx q[0];
rz(3.0359603) q[0];
rz(-pi) q[1];
rz(2.578031) q[2];
sx q[2];
rz(-2.5346018) q[2];
sx q[2];
rz(-2.3129472) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.48549451) q[1];
sx q[1];
rz(-1.6472677) q[1];
sx q[1];
rz(0.2904201) q[1];
rz(-pi) q[2];
rz(2.7518919) q[3];
sx q[3];
rz(-1.1003564) q[3];
sx q[3];
rz(-2.5999992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.303858) q[2];
sx q[2];
rz(-1.4759109) q[2];
sx q[2];
rz(-1.1676577) q[2];
rz(-1.6052823) q[3];
sx q[3];
rz(-1.4669908) q[3];
sx q[3];
rz(1.7355828) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7325571) q[0];
sx q[0];
rz(-1.5719825) q[0];
sx q[0];
rz(-2.4107966) q[0];
rz(2.2413975) q[1];
sx q[1];
rz(-2.3370445) q[1];
sx q[1];
rz(-2.3866167) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8687815) q[0];
sx q[0];
rz(-0.73505721) q[0];
sx q[0];
rz(-1.789577) q[0];
rz(-pi) q[1];
rz(0.99004284) q[2];
sx q[2];
rz(-0.92951894) q[2];
sx q[2];
rz(-2.5255447) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8440486) q[1];
sx q[1];
rz(-2.5522759) q[1];
sx q[1];
rz(-2.938016) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.8927535) q[3];
sx q[3];
rz(-2.4340981) q[3];
sx q[3];
rz(-1.0196109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3770611) q[2];
sx q[2];
rz(-1.3724644) q[2];
sx q[2];
rz(2.5047452) q[2];
rz(-2.8751255) q[3];
sx q[3];
rz(-2.0839432) q[3];
sx q[3];
rz(-1.586097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4144142) q[0];
sx q[0];
rz(-2.0120912) q[0];
sx q[0];
rz(1.138858) q[0];
rz(0.75421929) q[1];
sx q[1];
rz(-2.8051839) q[1];
sx q[1];
rz(-0.019502217) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0560023) q[0];
sx q[0];
rz(-1.3645571) q[0];
sx q[0];
rz(-3.0598559) q[0];
x q[1];
rz(-0.15140622) q[2];
sx q[2];
rz(-1.3726241) q[2];
sx q[2];
rz(2.1422447) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.86917415) q[1];
sx q[1];
rz(-0.65009102) q[1];
sx q[1];
rz(-1.8723349) q[1];
rz(-pi) q[2];
rz(2.9762514) q[3];
sx q[3];
rz(-1.2623566) q[3];
sx q[3];
rz(0.43499085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.1372244) q[2];
sx q[2];
rz(-1.4164111) q[2];
sx q[2];
rz(-0.84890378) q[2];
rz(-2.7539339) q[3];
sx q[3];
rz(-2.013423) q[3];
sx q[3];
rz(-1.5415812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7983109) q[0];
sx q[0];
rz(-0.16769519) q[0];
sx q[0];
rz(2.6570901) q[0];
rz(1.7548521) q[1];
sx q[1];
rz(-1.7157028) q[1];
sx q[1];
rz(-1.1482931) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0120221) q[0];
sx q[0];
rz(-2.9992933) q[0];
sx q[0];
rz(-2.9905031) q[0];
rz(0.30015595) q[2];
sx q[2];
rz(-0.65320063) q[2];
sx q[2];
rz(0.73210994) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.82603589) q[1];
sx q[1];
rz(-2.450374) q[1];
sx q[1];
rz(1.300699) q[1];
x q[2];
rz(-2.792312) q[3];
sx q[3];
rz(-1.619907) q[3];
sx q[3];
rz(1.9625361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.91223532) q[2];
sx q[2];
rz(-1.8476013) q[2];
sx q[2];
rz(-0.36515507) q[2];
rz(0.12864104) q[3];
sx q[3];
rz(-1.235685) q[3];
sx q[3];
rz(-2.685759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1098332) q[0];
sx q[0];
rz(-0.84072996) q[0];
sx q[0];
rz(1.6054556) q[0];
rz(-2.1784492) q[1];
sx q[1];
rz(-1.8704725) q[1];
sx q[1];
rz(2.0830547) q[1];
rz(-1.1744432) q[2];
sx q[2];
rz(-0.94559961) q[2];
sx q[2];
rz(-0.0018975817) q[2];
rz(-2.1309489) q[3];
sx q[3];
rz(-1.5394566) q[3];
sx q[3];
rz(1.4395366) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
