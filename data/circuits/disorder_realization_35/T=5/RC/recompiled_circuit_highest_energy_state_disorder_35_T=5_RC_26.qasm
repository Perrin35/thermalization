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
rz(-3.0776032) q[0];
sx q[0];
rz(-0.912323) q[0];
sx q[0];
rz(1.3349226) q[0];
rz(-1.1420684) q[1];
sx q[1];
rz(-0.51841441) q[1];
sx q[1];
rz(-1.4919182) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27150422) q[0];
sx q[0];
rz(-0.056500204) q[0];
sx q[0];
rz(-1.2369878) q[0];
x q[1];
rz(-0.9426192) q[2];
sx q[2];
rz(-2.16428) q[2];
sx q[2];
rz(-1.3243937) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3033437) q[1];
sx q[1];
rz(-2.456899) q[1];
sx q[1];
rz(2.0239023) q[1];
rz(-0.033871058) q[3];
sx q[3];
rz(-1.3283786) q[3];
sx q[3];
rz(-2.9592379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3492744) q[2];
sx q[2];
rz(-2.8742542) q[2];
sx q[2];
rz(3.086669) q[2];
rz(-1.6867636) q[3];
sx q[3];
rz(-0.39595404) q[3];
sx q[3];
rz(-1.6247862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9216264) q[0];
sx q[0];
rz(-1.3082137) q[0];
sx q[0];
rz(2.8935905) q[0];
rz(1.2451046) q[1];
sx q[1];
rz(-2.4047132) q[1];
sx q[1];
rz(-1.2287593) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4061444) q[0];
sx q[0];
rz(-1.7559253) q[0];
sx q[0];
rz(-1.6044751) q[0];
x q[1];
rz(2.6284559) q[2];
sx q[2];
rz(-2.0433132) q[2];
sx q[2];
rz(-3.0233011) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9039624) q[1];
sx q[1];
rz(-1.3591237) q[1];
sx q[1];
rz(-0.71118506) q[1];
rz(-pi) q[2];
rz(-0.13350962) q[3];
sx q[3];
rz(-0.4236246) q[3];
sx q[3];
rz(2.6299919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9241141) q[2];
sx q[2];
rz(-2.3096297) q[2];
sx q[2];
rz(0.58383101) q[2];
rz(0.42932388) q[3];
sx q[3];
rz(-1.9177633) q[3];
sx q[3];
rz(1.0113641) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9871224) q[0];
sx q[0];
rz(-2.7543289) q[0];
sx q[0];
rz(1.467147) q[0];
rz(-1.6636498) q[1];
sx q[1];
rz(-1.4091622) q[1];
sx q[1];
rz(2.0416226) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8407134) q[0];
sx q[0];
rz(-1.7689118) q[0];
sx q[0];
rz(1.677657) q[0];
x q[1];
rz(2.2713678) q[2];
sx q[2];
rz(-0.43284518) q[2];
sx q[2];
rz(-2.3926122) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8588455) q[1];
sx q[1];
rz(-1.9044151) q[1];
sx q[1];
rz(-2.7360914) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.93020029) q[3];
sx q[3];
rz(-0.67852321) q[3];
sx q[3];
rz(-2.2212706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.35525068) q[2];
sx q[2];
rz(-2.2189271) q[2];
sx q[2];
rz(-1.5938909) q[2];
rz(-1.330438) q[3];
sx q[3];
rz(-1.7682313) q[3];
sx q[3];
rz(1.0734585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7610157) q[0];
sx q[0];
rz(-1.8208193) q[0];
sx q[0];
rz(0.15783489) q[0];
rz(0.24179587) q[1];
sx q[1];
rz(-2.7561185) q[1];
sx q[1];
rz(-0.48113021) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24836981) q[0];
sx q[0];
rz(-2.4241894) q[0];
sx q[0];
rz(-3.0068842) q[0];
x q[1];
rz(-2.782576) q[2];
sx q[2];
rz(-1.8098157) q[2];
sx q[2];
rz(-2.0109107) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1222904) q[1];
sx q[1];
rz(-1.230352) q[1];
sx q[1];
rz(2.2346418) q[1];
rz(2.9322196) q[3];
sx q[3];
rz(-1.8518157) q[3];
sx q[3];
rz(-3.0176153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2699997) q[2];
sx q[2];
rz(-2.3303878) q[2];
sx q[2];
rz(-0.2743741) q[2];
rz(-1.4659878) q[3];
sx q[3];
rz(-1.2803187) q[3];
sx q[3];
rz(0.93332851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45193732) q[0];
sx q[0];
rz(-0.87757293) q[0];
sx q[0];
rz(2.080132) q[0];
rz(-2.6929216) q[1];
sx q[1];
rz(-2.6866388) q[1];
sx q[1];
rz(1.4400858) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3502304) q[0];
sx q[0];
rz(-1.4141948) q[0];
sx q[0];
rz(0.76275218) q[0];
rz(-1.1665197) q[2];
sx q[2];
rz(-1.6875522) q[2];
sx q[2];
rz(-1.5991885) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0764112) q[1];
sx q[1];
rz(-1.59756) q[1];
sx q[1];
rz(-1.8288495) q[1];
x q[2];
rz(2.1813356) q[3];
sx q[3];
rz(-1.4818076) q[3];
sx q[3];
rz(1.8833835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7777286) q[2];
sx q[2];
rz(-1.2960459) q[2];
sx q[2];
rz(-1.9484005) q[2];
rz(-2.7423972) q[3];
sx q[3];
rz(-1.4123071) q[3];
sx q[3];
rz(-3.1138368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4611918) q[0];
sx q[0];
rz(-1.2874648) q[0];
sx q[0];
rz(0.68049085) q[0];
rz(-0.83241278) q[1];
sx q[1];
rz(-1.2544371) q[1];
sx q[1];
rz(-0.15636538) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2966317) q[0];
sx q[0];
rz(-0.28066844) q[0];
sx q[0];
rz(0.190221) q[0];
rz(0.6758718) q[2];
sx q[2];
rz(-0.9524494) q[2];
sx q[2];
rz(-2.4175274) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.029441) q[1];
sx q[1];
rz(-1.2425155) q[1];
sx q[1];
rz(0.29354696) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0632319) q[3];
sx q[3];
rz(-2.8939191) q[3];
sx q[3];
rz(1.4192691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6303595) q[2];
sx q[2];
rz(-2.4801621) q[2];
sx q[2];
rz(-2.2516001) q[2];
rz(2.6651799) q[3];
sx q[3];
rz(-0.92372957) q[3];
sx q[3];
rz(0.43058968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-0.42590672) q[0];
sx q[0];
rz(-1.3193193) q[0];
sx q[0];
rz(2.529378) q[0];
rz(1.7828434) q[1];
sx q[1];
rz(-2.3814059) q[1];
sx q[1];
rz(-0.30119687) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8071547) q[0];
sx q[0];
rz(-0.2858735) q[0];
sx q[0];
rz(1.8854333) q[0];
rz(-pi) q[1];
rz(2.2280548) q[2];
sx q[2];
rz(-1.871167) q[2];
sx q[2];
rz(-2.0180839) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0665022) q[1];
sx q[1];
rz(-1.9478223) q[1];
sx q[1];
rz(2.6860363) q[1];
rz(-pi) q[2];
rz(2.878827) q[3];
sx q[3];
rz(-1.802236) q[3];
sx q[3];
rz(0.48514807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.23413868) q[2];
sx q[2];
rz(-0.80628866) q[2];
sx q[2];
rz(-1.0478728) q[2];
rz(-0.35510865) q[3];
sx q[3];
rz(-2.5706048) q[3];
sx q[3];
rz(-0.6855489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.909914) q[0];
sx q[0];
rz(-2.7993918) q[0];
sx q[0];
rz(-2.8177596) q[0];
rz(-2.553885) q[1];
sx q[1];
rz(-0.67981845) q[1];
sx q[1];
rz(-2.8388265) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2229109) q[0];
sx q[0];
rz(-0.41707539) q[0];
sx q[0];
rz(-2.4757451) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9043973) q[2];
sx q[2];
rz(-1.6182119) q[2];
sx q[2];
rz(0.42979017) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.98339048) q[1];
sx q[1];
rz(-1.2081971) q[1];
sx q[1];
rz(2.1561978) q[1];
rz(-pi) q[2];
rz(-2.8852194) q[3];
sx q[3];
rz(-1.1629761) q[3];
sx q[3];
rz(-2.2373509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4964464) q[2];
sx q[2];
rz(-1.030693) q[2];
sx q[2];
rz(-2.5808064) q[2];
rz(1.0049817) q[3];
sx q[3];
rz(-1.2882261) q[3];
sx q[3];
rz(1.0952449) q[3];
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
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1110558) q[0];
sx q[0];
rz(-0.66937864) q[0];
sx q[0];
rz(-2.3916767) q[0];
rz(2.7610682) q[1];
sx q[1];
rz(-0.98697248) q[1];
sx q[1];
rz(-0.97533018) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.153424) q[0];
sx q[0];
rz(-2.1200709) q[0];
sx q[0];
rz(2.2934521) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.777095) q[2];
sx q[2];
rz(-0.28536428) q[2];
sx q[2];
rz(-2.9138164) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6187894) q[1];
sx q[1];
rz(-2.0674043) q[1];
sx q[1];
rz(1.2148395) q[1];
x q[2];
rz(-2.1998243) q[3];
sx q[3];
rz(-2.7745651) q[3];
sx q[3];
rz(-2.0218771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9425977) q[2];
sx q[2];
rz(-2.2364605) q[2];
sx q[2];
rz(1.0814166) q[2];
rz(-0.069247581) q[3];
sx q[3];
rz(-1.705575) q[3];
sx q[3];
rz(0.92250219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13114318) q[0];
sx q[0];
rz(-0.32587019) q[0];
sx q[0];
rz(2.8358054) q[0];
rz(-0.30793134) q[1];
sx q[1];
rz(-1.6653857) q[1];
sx q[1];
rz(1.1526795) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9891805) q[0];
sx q[0];
rz(-2.3539641) q[0];
sx q[0];
rz(-2.7680725) q[0];
rz(-0.20736097) q[2];
sx q[2];
rz(-0.93610686) q[2];
sx q[2];
rz(2.7864252) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.88886966) q[1];
sx q[1];
rz(-0.53045853) q[1];
sx q[1];
rz(2.2525851) q[1];
rz(-pi) q[2];
rz(1.3373371) q[3];
sx q[3];
rz(-0.065107927) q[3];
sx q[3];
rz(-0.59681915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.61047381) q[2];
sx q[2];
rz(-1.0627154) q[2];
sx q[2];
rz(-1.4531685) q[2];
rz(0.48062634) q[3];
sx q[3];
rz(-0.56699816) q[3];
sx q[3];
rz(1.3948729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4881445) q[0];
sx q[0];
rz(-1.5339889) q[0];
sx q[0];
rz(1.4657159) q[0];
rz(-2.0987971) q[1];
sx q[1];
rz(-1.4029618) q[1];
sx q[1];
rz(-0.89697368) q[1];
rz(2.6063812) q[2];
sx q[2];
rz(-1.8028529) q[2];
sx q[2];
rz(0.48446083) q[2];
rz(0.35814169) q[3];
sx q[3];
rz(-1.5325357) q[3];
sx q[3];
rz(1.2383133) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
