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
rz(0.056869153) q[0];
sx q[0];
rz(2.9480204) q[0];
sx q[0];
rz(9.8326346) q[0];
rz(-0.017539311) q[1];
sx q[1];
rz(4.3932228) q[1];
sx q[1];
rz(7.8235758) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46288908) q[0];
sx q[0];
rz(-0.4400095) q[0];
sx q[0];
rz(1.4567039) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8694334) q[2];
sx q[2];
rz(-2.6725884) q[2];
sx q[2];
rz(2.9763165) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0527264) q[1];
sx q[1];
rz(-1.3590993) q[1];
sx q[1];
rz(2.746341) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2188721) q[3];
sx q[3];
rz(-2.2790103) q[3];
sx q[3];
rz(3.129209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.5876329) q[2];
sx q[2];
rz(-2.0822058) q[2];
sx q[2];
rz(2.9239192) q[2];
rz(-0.31146464) q[3];
sx q[3];
rz(-2.5296827) q[3];
sx q[3];
rz(-1.1658839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4199453) q[0];
sx q[0];
rz(-2.8657275) q[0];
sx q[0];
rz(-2.4496147) q[0];
rz(2.3852589) q[1];
sx q[1];
rz(-0.5223918) q[1];
sx q[1];
rz(-1.5501032) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4645351) q[0];
sx q[0];
rz(-1.4486396) q[0];
sx q[0];
rz(-2.4273901) q[0];
x q[1];
rz(1.9687551) q[2];
sx q[2];
rz(-1.6670456) q[2];
sx q[2];
rz(-0.41589662) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8263232) q[1];
sx q[1];
rz(-2.9401952) q[1];
sx q[1];
rz(-2.4884495) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0106929) q[3];
sx q[3];
rz(-2.0030795) q[3];
sx q[3];
rz(-3.0280857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.99593607) q[2];
sx q[2];
rz(-2.9684976) q[2];
sx q[2];
rz(-0.23925979) q[2];
rz(-1.650882) q[3];
sx q[3];
rz(-1.2942634) q[3];
sx q[3];
rz(-1.2283121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8656798) q[0];
sx q[0];
rz(-0.90621197) q[0];
sx q[0];
rz(3.1269585) q[0];
rz(-2.2485661) q[1];
sx q[1];
rz(-2.1664186) q[1];
sx q[1];
rz(2.9258974) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53447002) q[0];
sx q[0];
rz(-1.5661582) q[0];
sx q[0];
rz(1.6220553) q[0];
rz(2.7095238) q[2];
sx q[2];
rz(-2.4347907) q[2];
sx q[2];
rz(1.9343073) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.025192997) q[1];
sx q[1];
rz(-1.4433675) q[1];
sx q[1];
rz(-0.23372948) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.448378) q[3];
sx q[3];
rz(-2.1395464) q[3];
sx q[3];
rz(2.5162987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0821685) q[2];
sx q[2];
rz(-1.1801722) q[2];
sx q[2];
rz(0.99349418) q[2];
rz(-0.70837402) q[3];
sx q[3];
rz(-2.0230484) q[3];
sx q[3];
rz(0.12349252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74743903) q[0];
sx q[0];
rz(-0.49462947) q[0];
sx q[0];
rz(0.69394773) q[0];
rz(-2.8548062) q[1];
sx q[1];
rz(-1.6096121) q[1];
sx q[1];
rz(1.087711) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3932753) q[0];
sx q[0];
rz(-1.6946185) q[0];
sx q[0];
rz(1.0626777) q[0];
rz(-pi) q[1];
rz(0.7544341) q[2];
sx q[2];
rz(-2.4409119) q[2];
sx q[2];
rz(2.5452819) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.36632167) q[1];
sx q[1];
rz(-1.2676116) q[1];
sx q[1];
rz(2.1949185) q[1];
x q[2];
rz(-0.71432738) q[3];
sx q[3];
rz(-1.6052402) q[3];
sx q[3];
rz(-1.8689276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9168758) q[2];
sx q[2];
rz(-1.1563053) q[2];
sx q[2];
rz(-1.3678331) q[2];
rz(-1.2380838) q[3];
sx q[3];
rz(-1.5690469) q[3];
sx q[3];
rz(-2.9253173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5424159) q[0];
sx q[0];
rz(-1.2681862) q[0];
sx q[0];
rz(-2.5237778) q[0];
rz(-1.1082331) q[1];
sx q[1];
rz(-0.30312678) q[1];
sx q[1];
rz(-2.2584426) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9882433) q[0];
sx q[0];
rz(-0.81479615) q[0];
sx q[0];
rz(1.5809466) q[0];
x q[1];
rz(2.8034535) q[2];
sx q[2];
rz(-2.6827742) q[2];
sx q[2];
rz(2.3459319) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.42111629) q[1];
sx q[1];
rz(-1.5052374) q[1];
sx q[1];
rz(2.0620729) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0155025) q[3];
sx q[3];
rz(-1.2153373) q[3];
sx q[3];
rz(-2.773284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.21657476) q[2];
sx q[2];
rz(-2.8897132) q[2];
sx q[2];
rz(1.3612755) q[2];
rz(1.2733634) q[3];
sx q[3];
rz(-1.4342156) q[3];
sx q[3];
rz(0.98289615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.044416044) q[0];
sx q[0];
rz(-1.8955078) q[0];
sx q[0];
rz(2.521305) q[0];
rz(2.7445131) q[1];
sx q[1];
rz(-2.670791) q[1];
sx q[1];
rz(-1.9420067) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.202855) q[0];
sx q[0];
rz(-1.5917543) q[0];
sx q[0];
rz(-1.4644794) q[0];
x q[1];
rz(-1.3590089) q[2];
sx q[2];
rz(-1.4265991) q[2];
sx q[2];
rz(-1.5825001) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.37383693) q[1];
sx q[1];
rz(-1.7636313) q[1];
sx q[1];
rz(2.6111433) q[1];
x q[2];
rz(0.34440094) q[3];
sx q[3];
rz(-1.4592429) q[3];
sx q[3];
rz(1.2408181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8951796) q[2];
sx q[2];
rz(-1.2440888) q[2];
sx q[2];
rz(-2.4708774) q[2];
rz(-2.8504168) q[3];
sx q[3];
rz(-2.5305735) q[3];
sx q[3];
rz(0.62180716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42596844) q[0];
sx q[0];
rz(-1.427303) q[0];
sx q[0];
rz(-0.17844644) q[0];
rz(2.2715691) q[1];
sx q[1];
rz(-0.88302892) q[1];
sx q[1];
rz(-0.80284405) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74551302) q[0];
sx q[0];
rz(-0.95142309) q[0];
sx q[0];
rz(-1.5399163) q[0];
rz(0.3608256) q[2];
sx q[2];
rz(-0.61491167) q[2];
sx q[2];
rz(-0.71984839) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.72910129) q[1];
sx q[1];
rz(-1.8961398) q[1];
sx q[1];
rz(-0.71112432) q[1];
rz(0.25814806) q[3];
sx q[3];
rz(-0.95883639) q[3];
sx q[3];
rz(1.985509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.37658438) q[2];
sx q[2];
rz(-1.2976126) q[2];
sx q[2];
rz(2.2479748) q[2];
rz(-1.5417967) q[3];
sx q[3];
rz(-1.6518281) q[3];
sx q[3];
rz(2.5347575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97745085) q[0];
sx q[0];
rz(-0.41005382) q[0];
sx q[0];
rz(0.2151016) q[0];
rz(0.84456259) q[1];
sx q[1];
rz(-1.3203878) q[1];
sx q[1];
rz(0.79427687) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0117019) q[0];
sx q[0];
rz(-0.029057682) q[0];
sx q[0];
rz(-0.41883166) q[0];
rz(0.45659839) q[2];
sx q[2];
rz(-1.8482882) q[2];
sx q[2];
rz(-3.1268529) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2596742) q[1];
sx q[1];
rz(-1.6615903) q[1];
sx q[1];
rz(2.1114248) q[1];
x q[2];
rz(0.78171697) q[3];
sx q[3];
rz(-1.8060883) q[3];
sx q[3];
rz(2.5552251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8636785) q[2];
sx q[2];
rz(-1.2028376) q[2];
sx q[2];
rz(1.7318783) q[2];
rz(-2.388741) q[3];
sx q[3];
rz(-1.9949621) q[3];
sx q[3];
rz(-1.0021771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3442605) q[0];
sx q[0];
rz(-0.40818885) q[0];
sx q[0];
rz(0.97287384) q[0];
rz(1.6196039) q[1];
sx q[1];
rz(-1.7771959) q[1];
sx q[1];
rz(-2.5111228) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3909633) q[0];
sx q[0];
rz(-1.4328151) q[0];
sx q[0];
rz(2.9580381) q[0];
rz(-pi) q[1];
rz(3.1340808) q[2];
sx q[2];
rz(-1.4138616) q[2];
sx q[2];
rz(-0.73534009) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3296632) q[1];
sx q[1];
rz(-2.3721937) q[1];
sx q[1];
rz(0.071746221) q[1];
rz(-pi) q[2];
rz(-1.5077603) q[3];
sx q[3];
rz(-1.3975289) q[3];
sx q[3];
rz(1.2931371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.015532739) q[2];
sx q[2];
rz(-1.9893179) q[2];
sx q[2];
rz(-1.3014334) q[2];
rz(-1.5850916) q[3];
sx q[3];
rz(-0.82587487) q[3];
sx q[3];
rz(2.7263156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5179317) q[0];
sx q[0];
rz(-2.9908337) q[0];
sx q[0];
rz(2.6483722) q[0];
rz(-1.8863691) q[1];
sx q[1];
rz(-1.8938277) q[1];
sx q[1];
rz(-0.36466041) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55281524) q[0];
sx q[0];
rz(-1.7087666) q[0];
sx q[0];
rz(-1.9068043) q[0];
rz(-pi) q[1];
x q[1];
rz(0.45367806) q[2];
sx q[2];
rz(-0.96900207) q[2];
sx q[2];
rz(-1.2860677) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.2233878) q[1];
sx q[1];
rz(-2.3406696) q[1];
sx q[1];
rz(0.078322874) q[1];
rz(0.040382645) q[3];
sx q[3];
rz(-1.8869487) q[3];
sx q[3];
rz(-0.97163661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4517335) q[2];
sx q[2];
rz(-1.006459) q[2];
sx q[2];
rz(-0.9922007) q[2];
rz(2.774488) q[3];
sx q[3];
rz(-1.7016943) q[3];
sx q[3];
rz(-0.67830694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51919666) q[0];
sx q[0];
rz(-1.4764897) q[0];
sx q[0];
rz(1.3523703) q[0];
rz(-2.2182111) q[1];
sx q[1];
rz(-0.8538178) q[1];
sx q[1];
rz(-1.1710844) q[1];
rz(1.3707163) q[2];
sx q[2];
rz(-1.4309819) q[2];
sx q[2];
rz(-0.33071721) q[2];
rz(2.1711849) q[3];
sx q[3];
rz(-0.98894377) q[3];
sx q[3];
rz(0.43882216) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
