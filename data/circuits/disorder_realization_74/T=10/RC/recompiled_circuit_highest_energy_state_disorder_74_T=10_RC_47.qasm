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
rz(0.23298921) q[0];
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
rz(1.6132076) q[0];
sx q[0];
rz(-2.9068156) q[0];
sx q[0];
rz(1.6260207) q[0];
x q[1];
rz(-0.56143151) q[2];
sx q[2];
rz(-1.5806287) q[2];
sx q[2];
rz(-1.9468284) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0488551) q[1];
sx q[1];
rz(-1.6462781) q[1];
sx q[1];
rz(1.9239747) q[1];
x q[2];
rz(1.2337716) q[3];
sx q[3];
rz(-2.153018) q[3];
sx q[3];
rz(-1.9523417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.69922525) q[2];
sx q[2];
rz(-2.0825601) q[2];
sx q[2];
rz(3.0287058) q[2];
rz(2.8804307) q[3];
sx q[3];
rz(-1.7966725) q[3];
sx q[3];
rz(1.2507218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3498822) q[0];
sx q[0];
rz(-0.078502027) q[0];
sx q[0];
rz(3.0872524) q[0];
rz(0.21866523) q[1];
sx q[1];
rz(-1.4726787) q[1];
sx q[1];
rz(-2.7770619) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.88663) q[0];
sx q[0];
rz(-0.97575649) q[0];
sx q[0];
rz(-0.9592077) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.40016115) q[2];
sx q[2];
rz(-2.1605943) q[2];
sx q[2];
rz(1.9831744) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5082701) q[1];
sx q[1];
rz(-2.5145217) q[1];
sx q[1];
rz(-0.12992628) q[1];
x q[2];
rz(0.18192911) q[3];
sx q[3];
rz(-1.9046138) q[3];
sx q[3];
rz(-0.16063375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9008122) q[2];
sx q[2];
rz(-1.6996926) q[2];
sx q[2];
rz(2.0330632) q[2];
rz(0.19567868) q[3];
sx q[3];
rz(-1.9340197) q[3];
sx q[3];
rz(1.4669363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(2.7473258) q[0];
sx q[0];
rz(-2.1013923) q[0];
sx q[0];
rz(2.105383) q[0];
rz(0.58468435) q[1];
sx q[1];
rz(-1.6744637) q[1];
sx q[1];
rz(2.5856957) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3506539) q[0];
sx q[0];
rz(-2.2430111) q[0];
sx q[0];
rz(-2.1124798) q[0];
x q[1];
rz(2.0228407) q[2];
sx q[2];
rz(-3.0556745) q[2];
sx q[2];
rz(-2.265386) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8396254) q[1];
sx q[1];
rz(-0.71528597) q[1];
sx q[1];
rz(1.7710745) q[1];
rz(1.317765) q[3];
sx q[3];
rz(-2.0076499) q[3];
sx q[3];
rz(1.4405314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3802152) q[2];
sx q[2];
rz(-2.9351202) q[2];
sx q[2];
rz(-1.9492487) q[2];
rz(-3.0377667) q[3];
sx q[3];
rz(-1.9793648) q[3];
sx q[3];
rz(1.9942358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.035606774) q[0];
sx q[0];
rz(-0.62653956) q[0];
sx q[0];
rz(-1.71738) q[0];
rz(2.4183938) q[1];
sx q[1];
rz(-0.23591787) q[1];
sx q[1];
rz(0.22044388) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28726668) q[0];
sx q[0];
rz(-1.1136855) q[0];
sx q[0];
rz(2.7316389) q[0];
x q[1];
rz(1.0102369) q[2];
sx q[2];
rz(-1.3187485) q[2];
sx q[2];
rz(-0.28365669) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6666777) q[1];
sx q[1];
rz(-1.6834458) q[1];
sx q[1];
rz(-1.2126318) q[1];
rz(-pi) q[2];
rz(-1.1573767) q[3];
sx q[3];
rz(-1.2261021) q[3];
sx q[3];
rz(2.7935296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7580938) q[2];
sx q[2];
rz(-1.8269962) q[2];
sx q[2];
rz(-0.51551762) q[2];
rz(3.0564485) q[3];
sx q[3];
rz(-0.65535039) q[3];
sx q[3];
rz(2.3084124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4715217) q[0];
sx q[0];
rz(-2.9586198) q[0];
sx q[0];
rz(-0.66437379) q[0];
rz(-0.5270671) q[1];
sx q[1];
rz(-1.138569) q[1];
sx q[1];
rz(1.1767496) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7648286) q[0];
sx q[0];
rz(-1.6559147) q[0];
sx q[0];
rz(3.1133988) q[0];
rz(-2.2559153) q[2];
sx q[2];
rz(-1.9373477) q[2];
sx q[2];
rz(-1.4057856) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6709395) q[1];
sx q[1];
rz(-1.8090418) q[1];
sx q[1];
rz(-1.5508503) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.845729) q[3];
sx q[3];
rz(-2.5524013) q[3];
sx q[3];
rz(-1.2146035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.027017) q[2];
sx q[2];
rz(-2.9288374) q[2];
sx q[2];
rz(-2.2034755) q[2];
rz(1.045687) q[3];
sx q[3];
rz(-0.18200471) q[3];
sx q[3];
rz(-0.81030455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(2.3351347) q[0];
sx q[0];
rz(-1.198575) q[0];
sx q[0];
rz(-1.2499811) q[0];
rz(-0.20420034) q[1];
sx q[1];
rz(-0.67664346) q[1];
sx q[1];
rz(2.2883889) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1565022) q[0];
sx q[0];
rz(-1.2285875) q[0];
sx q[0];
rz(1.8746822) q[0];
rz(-2.5218042) q[2];
sx q[2];
rz(-1.0519805) q[2];
sx q[2];
rz(2.3960115) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.051999747) q[1];
sx q[1];
rz(-1.9885364) q[1];
sx q[1];
rz(-2.3485648) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.12746396) q[3];
sx q[3];
rz(-1.9098305) q[3];
sx q[3];
rz(-1.6515428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0975254) q[2];
sx q[2];
rz(-1.3419469) q[2];
sx q[2];
rz(-1.1996783) q[2];
rz(-2.1357644) q[3];
sx q[3];
rz(-0.89322105) q[3];
sx q[3];
rz(-2.306166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47578874) q[0];
sx q[0];
rz(-2.8732193) q[0];
sx q[0];
rz(2.1589808) q[0];
rz(1.8334552) q[1];
sx q[1];
rz(-1.2160701) q[1];
sx q[1];
rz(1.4391724) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60123283) q[0];
sx q[0];
rz(-1.2315016) q[0];
sx q[0];
rz(-1.220006) q[0];
x q[1];
rz(0.077567049) q[2];
sx q[2];
rz(-1.0500538) q[2];
sx q[2];
rz(-2.5375053) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.01165302) q[1];
sx q[1];
rz(-2.2352152) q[1];
sx q[1];
rz(2.7205977) q[1];
rz(-pi) q[2];
rz(-1.7129219) q[3];
sx q[3];
rz(-1.097659) q[3];
sx q[3];
rz(0.17156916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9787489) q[2];
sx q[2];
rz(-0.89357251) q[2];
sx q[2];
rz(-0.54235512) q[2];
rz(0.98313037) q[3];
sx q[3];
rz(-1.1636846) q[3];
sx q[3];
rz(0.94016176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5777609) q[0];
sx q[0];
rz(-1.4262119) q[0];
sx q[0];
rz(-2.3610709) q[0];
rz(1.1753987) q[1];
sx q[1];
rz(-0.54882097) q[1];
sx q[1];
rz(-1.0343879) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9517453) q[0];
sx q[0];
rz(-1.4035514) q[0];
sx q[0];
rz(2.7664004) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.081806094) q[2];
sx q[2];
rz(-1.713364) q[2];
sx q[2];
rz(0.35065251) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8966297) q[1];
sx q[1];
rz(-1.9707929) q[1];
sx q[1];
rz(-1.4463615) q[1];
rz(-pi) q[2];
x q[2];
rz(0.54776056) q[3];
sx q[3];
rz(-1.3410853) q[3];
sx q[3];
rz(-1.3198157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.15022755) q[2];
sx q[2];
rz(-1.4069724) q[2];
sx q[2];
rz(0.11037174) q[2];
rz(-1.2982093) q[3];
sx q[3];
rz(-1.7896264) q[3];
sx q[3];
rz(-2.8492294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89123911) q[0];
sx q[0];
rz(-1.013101) q[0];
sx q[0];
rz(-1.2497586) q[0];
rz(0.091014422) q[1];
sx q[1];
rz(-2.4200771) q[1];
sx q[1];
rz(0.7116085) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7994493) q[0];
sx q[0];
rz(-0.79252386) q[0];
sx q[0];
rz(-1.5802797) q[0];
rz(-pi) q[1];
rz(-0.83699147) q[2];
sx q[2];
rz(-2.2185433) q[2];
sx q[2];
rz(-0.68488065) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8636057) q[1];
sx q[1];
rz(-0.58574876) q[1];
sx q[1];
rz(0.65237712) q[1];
rz(-pi) q[2];
rz(-1.9125838) q[3];
sx q[3];
rz(-2.2967489) q[3];
sx q[3];
rz(2.1330698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8083501) q[2];
sx q[2];
rz(-1.658172) q[2];
sx q[2];
rz(-1.1727772) q[2];
rz(1.6138389) q[3];
sx q[3];
rz(-2.0634191) q[3];
sx q[3];
rz(-3.0465904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.461819) q[0];
sx q[0];
rz(-3.0175896) q[0];
sx q[0];
rz(1.5661731) q[0];
rz(0.35789403) q[1];
sx q[1];
rz(-1.0708258) q[1];
sx q[1];
rz(-0.90248743) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76644635) q[0];
sx q[0];
rz(-1.8286819) q[0];
sx q[0];
rz(-0.93574406) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4975791) q[2];
sx q[2];
rz(-1.5884382) q[2];
sx q[2];
rz(1.0746852) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.27943006) q[1];
sx q[1];
rz(-1.2365274) q[1];
sx q[1];
rz(-2.3390655) q[1];
rz(-pi) q[2];
x q[2];
rz(0.087335056) q[3];
sx q[3];
rz(-2.119434) q[3];
sx q[3];
rz(-0.61591301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0365399) q[2];
sx q[2];
rz(-2.5532711) q[2];
sx q[2];
rz(-1.6416637) q[2];
rz(-0.52122742) q[3];
sx q[3];
rz(-0.14385496) q[3];
sx q[3];
rz(2.1874793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4317516) q[0];
sx q[0];
rz(-0.7702282) q[0];
sx q[0];
rz(-1.7256398) q[0];
rz(-0.00090986666) q[1];
sx q[1];
rz(-1.6699427) q[1];
sx q[1];
rz(-1.4572399) q[1];
rz(0.098401423) q[2];
sx q[2];
rz(-2.142061) q[2];
sx q[2];
rz(-0.67666035) q[2];
rz(-0.74827452) q[3];
sx q[3];
rz(-2.1088441) q[3];
sx q[3];
rz(0.90499457) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
