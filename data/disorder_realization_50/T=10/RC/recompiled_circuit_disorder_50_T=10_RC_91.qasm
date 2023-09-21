OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.87941909) q[0];
sx q[0];
rz(4.8882422) q[0];
sx q[0];
rz(12.56765) q[0];
rz(1.4446422) q[1];
sx q[1];
rz(-1.1029707) q[1];
sx q[1];
rz(-0.77494088) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2389195) q[0];
sx q[0];
rz(-2.9905149) q[0];
sx q[0];
rz(2.8047049) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.84471976) q[2];
sx q[2];
rz(-1.2436927) q[2];
sx q[2];
rz(-2.6127882) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.85841132) q[1];
sx q[1];
rz(-1.6106669) q[1];
sx q[1];
rz(-1.4285425) q[1];
rz(-pi) q[2];
rz(-1.418387) q[3];
sx q[3];
rz(-1.2890352) q[3];
sx q[3];
rz(-1.2581274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0455735) q[2];
sx q[2];
rz(-2.0298268) q[2];
sx q[2];
rz(1.1958896) q[2];
rz(-1.1536417) q[3];
sx q[3];
rz(-0.78912815) q[3];
sx q[3];
rz(1.7529863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.4202704) q[0];
sx q[0];
rz(-0.35636154) q[0];
sx q[0];
rz(-1.2715682) q[0];
rz(1.0999854) q[1];
sx q[1];
rz(-2.0309235) q[1];
sx q[1];
rz(1.3756479) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.401706) q[0];
sx q[0];
rz(-1.4835843) q[0];
sx q[0];
rz(0.4102104) q[0];
x q[1];
rz(3.0963418) q[2];
sx q[2];
rz(-1.3327193) q[2];
sx q[2];
rz(2.0766052) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.34054204) q[1];
sx q[1];
rz(-1.935563) q[1];
sx q[1];
rz(1.7629452) q[1];
rz(-pi) q[2];
x q[2];
rz(0.35956412) q[3];
sx q[3];
rz(-1.0457977) q[3];
sx q[3];
rz(-2.0294702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9282844) q[2];
sx q[2];
rz(-1.3457315) q[2];
sx q[2];
rz(2.6518872) q[2];
rz(1.1335763) q[3];
sx q[3];
rz(-2.9462892) q[3];
sx q[3];
rz(1.0666696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60004822) q[0];
sx q[0];
rz(-1.8627889) q[0];
sx q[0];
rz(-0.24060732) q[0];
rz(-0.34257564) q[1];
sx q[1];
rz(-0.97476417) q[1];
sx q[1];
rz(-1.906357) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2152527) q[0];
sx q[0];
rz(-2.140576) q[0];
sx q[0];
rz(1.1220054) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7690976) q[2];
sx q[2];
rz(-1.3349581) q[2];
sx q[2];
rz(1.3553938) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8894549) q[1];
sx q[1];
rz(-2.7733907) q[1];
sx q[1];
rz(0.77106573) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5123185) q[3];
sx q[3];
rz(-0.73262239) q[3];
sx q[3];
rz(2.7936558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3774595) q[2];
sx q[2];
rz(-1.5185792) q[2];
sx q[2];
rz(-1.7049449) q[2];
rz(1.7403587) q[3];
sx q[3];
rz(-1.2652206) q[3];
sx q[3];
rz(-2.5333372) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.075994611) q[0];
sx q[0];
rz(-1.6125212) q[0];
sx q[0];
rz(-2.3377989) q[0];
rz(-0.94961387) q[1];
sx q[1];
rz(-1.4664374) q[1];
sx q[1];
rz(0.11985699) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0940995) q[0];
sx q[0];
rz(-2.1764628) q[0];
sx q[0];
rz(0.043461965) q[0];
rz(-1.4255964) q[2];
sx q[2];
rz(-1.2135266) q[2];
sx q[2];
rz(1.4841339) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4704628) q[1];
sx q[1];
rz(-1.5884591) q[1];
sx q[1];
rz(-1.8030241) q[1];
x q[2];
rz(-1.0355789) q[3];
sx q[3];
rz(-2.1018627) q[3];
sx q[3];
rz(2.5553184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.63311657) q[2];
sx q[2];
rz(-1.8957596) q[2];
sx q[2];
rz(-0.072337739) q[2];
rz(-0.37483254) q[3];
sx q[3];
rz(-0.65632498) q[3];
sx q[3];
rz(-1.9434631) q[3];
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
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4478093) q[0];
sx q[0];
rz(-0.57149514) q[0];
sx q[0];
rz(2.6547292) q[0];
rz(0.72987366) q[1];
sx q[1];
rz(-2.2327773) q[1];
sx q[1];
rz(-1.9015076) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64802058) q[0];
sx q[0];
rz(-0.23952142) q[0];
sx q[0];
rz(-3.0859205) q[0];
x q[1];
rz(0.24057062) q[2];
sx q[2];
rz(-1.3153937) q[2];
sx q[2];
rz(-2.0672928) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6606635) q[1];
sx q[1];
rz(-0.96748057) q[1];
sx q[1];
rz(-1.6515345) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1754856) q[3];
sx q[3];
rz(-1.9786414) q[3];
sx q[3];
rz(0.84433489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1029677) q[2];
sx q[2];
rz(-0.90927783) q[2];
sx q[2];
rz(2.4385578) q[2];
rz(-1.7317584) q[3];
sx q[3];
rz(-1.3491646) q[3];
sx q[3];
rz(2.9838802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(-0.96419656) q[0];
sx q[0];
rz(-1.8494158) q[0];
sx q[0];
rz(0.99037209) q[0];
rz(0.05274996) q[1];
sx q[1];
rz(-0.92664781) q[1];
sx q[1];
rz(1.4809158) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4303362) q[0];
sx q[0];
rz(-2.7529786) q[0];
sx q[0];
rz(-1.7646164) q[0];
rz(1.5724206) q[2];
sx q[2];
rz(-2.6964028) q[2];
sx q[2];
rz(-2.0507398) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8112091) q[1];
sx q[1];
rz(-0.57302176) q[1];
sx q[1];
rz(1.9622383) q[1];
rz(-pi) q[2];
rz(-0.048599676) q[3];
sx q[3];
rz(-0.89266333) q[3];
sx q[3];
rz(-2.5130659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.0084373077) q[2];
sx q[2];
rz(-1.6494273) q[2];
sx q[2];
rz(-0.56224242) q[2];
rz(1.0605313) q[3];
sx q[3];
rz(-2.3980467) q[3];
sx q[3];
rz(2.8806768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9437207) q[0];
sx q[0];
rz(-1.3523538) q[0];
sx q[0];
rz(-1.4690171) q[0];
rz(-2.127227) q[1];
sx q[1];
rz(-2.1140153) q[1];
sx q[1];
rz(-1.8168824) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0425456) q[0];
sx q[0];
rz(-1.5481661) q[0];
sx q[0];
rz(-3.0142473) q[0];
rz(-2.1340738) q[2];
sx q[2];
rz(-1.129732) q[2];
sx q[2];
rz(1.2457459) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.99609375) q[1];
sx q[1];
rz(-0.56963527) q[1];
sx q[1];
rz(2.3949404) q[1];
rz(-pi) q[2];
rz(2.0539829) q[3];
sx q[3];
rz(-2.4251068) q[3];
sx q[3];
rz(-2.0487752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5614732) q[2];
sx q[2];
rz(-1.7529528) q[2];
sx q[2];
rz(0.44357792) q[2];
rz(0.94868547) q[3];
sx q[3];
rz(-1.1921927) q[3];
sx q[3];
rz(2.366812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.401944) q[0];
sx q[0];
rz(-1.0111324) q[0];
sx q[0];
rz(-0.46052128) q[0];
rz(3.0415688) q[1];
sx q[1];
rz(-2.1477551) q[1];
sx q[1];
rz(1.2896279) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3279009) q[0];
sx q[0];
rz(-1.1192338) q[0];
sx q[0];
rz(2.5149462) q[0];
rz(-pi) q[1];
rz(-0.71596594) q[2];
sx q[2];
rz(-2.1858474) q[2];
sx q[2];
rz(1.3032608) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.928927) q[1];
sx q[1];
rz(-1.7576828) q[1];
sx q[1];
rz(-3.0287663) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7361705) q[3];
sx q[3];
rz(-2.227042) q[3];
sx q[3];
rz(-1.0429494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.148968) q[2];
sx q[2];
rz(-2.830539) q[2];
sx q[2];
rz(2.0588493) q[2];
rz(3.0454214) q[3];
sx q[3];
rz(-1.7008737) q[3];
sx q[3];
rz(1.2307897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37725317) q[0];
sx q[0];
rz(-2.9057673) q[0];
sx q[0];
rz(-2.8073231) q[0];
rz(1.224068) q[1];
sx q[1];
rz(-1.5609488) q[1];
sx q[1];
rz(-0.28265488) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44336549) q[0];
sx q[0];
rz(-2.3400314) q[0];
sx q[0];
rz(-2.482224) q[0];
x q[1];
rz(0.56844175) q[2];
sx q[2];
rz(-2.754515) q[2];
sx q[2];
rz(0.27075152) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0375367) q[1];
sx q[1];
rz(-1.1661342) q[1];
sx q[1];
rz(-1.3794823) q[1];
rz(2.2563124) q[3];
sx q[3];
rz(-1.4797398) q[3];
sx q[3];
rz(1.0202927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.21215542) q[2];
sx q[2];
rz(-1.9781457) q[2];
sx q[2];
rz(-0.7412509) q[2];
rz(0.50179982) q[3];
sx q[3];
rz(-1.8024249) q[3];
sx q[3];
rz(-2.5575976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(2.042645) q[0];
sx q[0];
rz(-0.89530033) q[0];
sx q[0];
rz(0.50869554) q[0];
rz(0.11518654) q[1];
sx q[1];
rz(-2.4337264) q[1];
sx q[1];
rz(2.4597816) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54036056) q[0];
sx q[0];
rz(-1.0245748) q[0];
sx q[0];
rz(2.9451314) q[0];
rz(-pi) q[1];
rz(2.4585312) q[2];
sx q[2];
rz(-2.1944754) q[2];
sx q[2];
rz(-1.1766441) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0444813) q[1];
sx q[1];
rz(-1.0978571) q[1];
sx q[1];
rz(0.60826917) q[1];
rz(-pi) q[2];
x q[2];
rz(0.82716771) q[3];
sx q[3];
rz(-0.70561545) q[3];
sx q[3];
rz(3.0527715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8455785) q[2];
sx q[2];
rz(-1.5073551) q[2];
sx q[2];
rz(0.34109035) q[2];
rz(-1.0836481) q[3];
sx q[3];
rz(-2.3458979) q[3];
sx q[3];
rz(-0.17102374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1682128) q[0];
sx q[0];
rz(-1.4773049) q[0];
sx q[0];
rz(2.1422577) q[0];
rz(-2.534261) q[1];
sx q[1];
rz(-0.9098396) q[1];
sx q[1];
rz(0.20914016) q[1];
rz(-2.091999) q[2];
sx q[2];
rz(-0.97432077) q[2];
sx q[2];
rz(-1.5885098) q[2];
rz(-0.46799216) q[3];
sx q[3];
rz(-2.7157126) q[3];
sx q[3];
rz(-1.4117905) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
