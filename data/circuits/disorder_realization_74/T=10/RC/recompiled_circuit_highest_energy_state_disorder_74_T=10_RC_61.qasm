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
rz(1.9407152) q[0];
sx q[0];
rz(-2.241029) q[0];
sx q[0];
rz(-2.9086034) q[0];
rz(1.6147344) q[1];
sx q[1];
rz(-1.072071) q[1];
sx q[1];
rz(2.022068) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6699864) q[0];
sx q[0];
rz(-1.8052088) q[0];
sx q[0];
rz(-0.013201518) q[0];
rz(-0.56143151) q[2];
sx q[2];
rz(-1.5806287) q[2];
sx q[2];
rz(1.1947643) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0488551) q[1];
sx q[1];
rz(-1.6462781) q[1];
sx q[1];
rz(1.2176179) q[1];
rz(-1.9078211) q[3];
sx q[3];
rz(-0.98857461) q[3];
sx q[3];
rz(1.9523417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.69922525) q[2];
sx q[2];
rz(-2.0825601) q[2];
sx q[2];
rz(3.0287058) q[2];
rz(-2.8804307) q[3];
sx q[3];
rz(-1.3449202) q[3];
sx q[3];
rz(1.2507218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3498822) q[0];
sx q[0];
rz(-0.078502027) q[0];
sx q[0];
rz(-3.0872524) q[0];
rz(0.21866523) q[1];
sx q[1];
rz(-1.668914) q[1];
sx q[1];
rz(-0.36453077) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69037777) q[0];
sx q[0];
rz(-1.0752819) q[0];
sx q[0];
rz(2.4507603) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1991992) q[2];
sx q[2];
rz(-1.2411054) q[2];
sx q[2];
rz(0.18131944) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6682525) q[1];
sx q[1];
rz(-2.1917731) q[1];
sx q[1];
rz(-1.6644068) q[1];
rz(-pi) q[2];
rz(-1.909799) q[3];
sx q[3];
rz(-1.7425797) q[3];
sx q[3];
rz(1.4703657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.24078044) q[2];
sx q[2];
rz(-1.6996926) q[2];
sx q[2];
rz(-2.0330632) q[2];
rz(0.19567868) q[3];
sx q[3];
rz(-1.9340197) q[3];
sx q[3];
rz(-1.6746563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3942669) q[0];
sx q[0];
rz(-2.1013923) q[0];
sx q[0];
rz(-1.0362097) q[0];
rz(0.58468435) q[1];
sx q[1];
rz(-1.4671289) q[1];
sx q[1];
rz(0.55589693) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42134684) q[0];
sx q[0];
rz(-1.9860391) q[0];
sx q[0];
rz(-0.74851997) q[0];
x q[1];
rz(2.0228407) q[2];
sx q[2];
rz(-0.085918203) q[2];
sx q[2];
rz(-0.8762067) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8396254) q[1];
sx q[1];
rz(-0.71528597) q[1];
sx q[1];
rz(-1.3705181) q[1];
x q[2];
rz(-0.44938748) q[3];
sx q[3];
rz(-1.7996178) q[3];
sx q[3];
rz(0.23923161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7613775) q[2];
sx q[2];
rz(-0.20647241) q[2];
sx q[2];
rz(1.192344) q[2];
rz(-0.10382593) q[3];
sx q[3];
rz(-1.9793648) q[3];
sx q[3];
rz(-1.9942358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.035606774) q[0];
sx q[0];
rz(-2.5150531) q[0];
sx q[0];
rz(-1.4242127) q[0];
rz(0.72319889) q[1];
sx q[1];
rz(-0.23591787) q[1];
sx q[1];
rz(2.9211488) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.854326) q[0];
sx q[0];
rz(-1.1136855) q[0];
sx q[0];
rz(2.7316389) q[0];
x q[1];
rz(2.0218684) q[2];
sx q[2];
rz(-2.5325518) q[2];
sx q[2];
rz(0.90897564) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6666777) q[1];
sx q[1];
rz(-1.6834458) q[1];
sx q[1];
rz(-1.9289609) q[1];
rz(-pi) q[2];
rz(1.9842159) q[3];
sx q[3];
rz(-1.2261021) q[3];
sx q[3];
rz(2.7935296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.38349884) q[2];
sx q[2];
rz(-1.8269962) q[2];
sx q[2];
rz(2.626075) q[2];
rz(-0.085144194) q[3];
sx q[3];
rz(-0.65535039) q[3];
sx q[3];
rz(-0.83318025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67007095) q[0];
sx q[0];
rz(-0.18297289) q[0];
sx q[0];
rz(-2.4772189) q[0];
rz(-0.5270671) q[1];
sx q[1];
rz(-2.0030237) q[1];
sx q[1];
rz(1.9648431) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3767641) q[0];
sx q[0];
rz(-1.6559147) q[0];
sx q[0];
rz(3.1133988) q[0];
rz(-pi) q[1];
rz(1.0254622) q[2];
sx q[2];
rz(-2.3787913) q[2];
sx q[2];
rz(-2.5632586) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5866111) q[1];
sx q[1];
rz(-0.23906318) q[1];
sx q[1];
rz(-0.081940941) q[1];
rz(-pi) q[2];
rz(1.2958636) q[3];
sx q[3];
rz(-2.5524013) q[3];
sx q[3];
rz(1.9269892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1145757) q[2];
sx q[2];
rz(-2.9288374) q[2];
sx q[2];
rz(2.2034755) q[2];
rz(1.045687) q[3];
sx q[3];
rz(-2.9595879) q[3];
sx q[3];
rz(0.81030455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3351347) q[0];
sx q[0];
rz(-1.9430176) q[0];
sx q[0];
rz(1.8916116) q[0];
rz(2.9373923) q[1];
sx q[1];
rz(-0.67664346) q[1];
sx q[1];
rz(2.2883889) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7366587) q[0];
sx q[0];
rz(-0.45365004) q[0];
sx q[0];
rz(2.4429295) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5218042) q[2];
sx q[2];
rz(-2.0896122) q[2];
sx q[2];
rz(2.3960115) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.051999747) q[1];
sx q[1];
rz(-1.9885364) q[1];
sx q[1];
rz(-0.79302782) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2248269) q[3];
sx q[3];
rz(-0.36133535) q[3];
sx q[3];
rz(-1.122235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0975254) q[2];
sx q[2];
rz(-1.7996457) q[2];
sx q[2];
rz(1.9419144) q[2];
rz(2.1357644) q[3];
sx q[3];
rz(-2.2483716) q[3];
sx q[3];
rz(0.83542663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6658039) q[0];
sx q[0];
rz(-2.8732193) q[0];
sx q[0];
rz(0.98261181) q[0];
rz(-1.3081374) q[1];
sx q[1];
rz(-1.9255226) q[1];
sx q[1];
rz(1.7024202) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60123283) q[0];
sx q[0];
rz(-1.2315016) q[0];
sx q[0];
rz(-1.220006) q[0];
x q[1];
rz(-3.0640256) q[2];
sx q[2];
rz(-2.0915389) q[2];
sx q[2];
rz(2.5375053) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.1299396) q[1];
sx q[1];
rz(-2.2352152) q[1];
sx q[1];
rz(-0.42099492) q[1];
rz(-1.4286707) q[3];
sx q[3];
rz(-1.097659) q[3];
sx q[3];
rz(2.9700235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9787489) q[2];
sx q[2];
rz(-0.89357251) q[2];
sx q[2];
rz(2.5992375) q[2];
rz(-0.98313037) q[3];
sx q[3];
rz(-1.977908) q[3];
sx q[3];
rz(-2.2014309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56383175) q[0];
sx q[0];
rz(-1.7153808) q[0];
sx q[0];
rz(2.3610709) q[0];
rz(-1.1753987) q[1];
sx q[1];
rz(-2.5927717) q[1];
sx q[1];
rz(-1.0343879) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18984737) q[0];
sx q[0];
rz(-1.4035514) q[0];
sx q[0];
rz(-2.7664004) q[0];
rz(-1.7138359) q[2];
sx q[2];
rz(-1.6517706) q[2];
sx q[2];
rz(1.2317927) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7670892) q[1];
sx q[1];
rz(-1.6853635) q[1];
sx q[1];
rz(2.7388045) q[1];
x q[2];
rz(-1.8381528) q[3];
sx q[3];
rz(-2.1026097) q[3];
sx q[3];
rz(-0.11296266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.15022755) q[2];
sx q[2];
rz(-1.4069724) q[2];
sx q[2];
rz(-3.0312209) q[2];
rz(-1.8433833) q[3];
sx q[3];
rz(-1.3519663) q[3];
sx q[3];
rz(0.29236326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89123911) q[0];
sx q[0];
rz(-1.013101) q[0];
sx q[0];
rz(-1.891834) q[0];
rz(0.091014422) q[1];
sx q[1];
rz(-0.72151557) q[1];
sx q[1];
rz(-0.7116085) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7994493) q[0];
sx q[0];
rz(-0.79252386) q[0];
sx q[0];
rz(-1.5802797) q[0];
rz(-pi) q[1];
rz(0.83699147) q[2];
sx q[2];
rz(-2.2185433) q[2];
sx q[2];
rz(0.68488065) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8636057) q[1];
sx q[1];
rz(-0.58574876) q[1];
sx q[1];
rz(0.65237712) q[1];
rz(1.2290088) q[3];
sx q[3];
rz(-0.84484378) q[3];
sx q[3];
rz(-2.1330698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.3332425) q[2];
sx q[2];
rz(-1.4834206) q[2];
sx q[2];
rz(-1.9688155) q[2];
rz(-1.6138389) q[3];
sx q[3];
rz(-1.0781735) q[3];
sx q[3];
rz(0.095002256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67977366) q[0];
sx q[0];
rz(-0.12400308) q[0];
sx q[0];
rz(-1.5661731) q[0];
rz(0.35789403) q[1];
sx q[1];
rz(-1.0708258) q[1];
sx q[1];
rz(-0.90248743) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99011496) q[0];
sx q[0];
rz(-0.95989908) q[0];
sx q[0];
rz(2.8249802) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.64401354) q[2];
sx q[2];
rz(-1.5884382) q[2];
sx q[2];
rz(-1.0746852) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1572307) q[1];
sx q[1];
rz(-2.2868726) q[1];
sx q[1];
rz(-2.6916719) q[1];
rz(-pi) q[2];
rz(1.020458) q[3];
sx q[3];
rz(-1.4963048) q[3];
sx q[3];
rz(2.1410774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1050528) q[2];
sx q[2];
rz(-0.5883216) q[2];
sx q[2];
rz(-1.6416637) q[2];
rz(-2.6203652) q[3];
sx q[3];
rz(-0.14385496) q[3];
sx q[3];
rz(0.95411333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
rz(0.70984107) q[0];
sx q[0];
rz(-2.3713645) q[0];
sx q[0];
rz(1.4159528) q[0];
rz(0.00090986666) q[1];
sx q[1];
rz(-1.4716499) q[1];
sx q[1];
rz(1.6843527) q[1];
rz(-1.7224689) q[2];
sx q[2];
rz(-2.5628452) q[2];
sx q[2];
rz(2.645523) q[2];
rz(-0.88738995) q[3];
sx q[3];
rz(-2.1947464) q[3];
sx q[3];
rz(-1.1098679) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
