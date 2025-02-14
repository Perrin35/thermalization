OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.5105628) q[0];
sx q[0];
rz(-0.50463843) q[0];
sx q[0];
rz(2.7161427) q[0];
rz(0.57234859) q[1];
sx q[1];
rz(-2.0444137) q[1];
sx q[1];
rz(-1.9414577) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.67681) q[0];
sx q[0];
rz(-1.3199782) q[0];
sx q[0];
rz(-1.5607911) q[0];
rz(2.8699371) q[2];
sx q[2];
rz(-1.9052398) q[2];
sx q[2];
rz(1.7890499) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8737826) q[1];
sx q[1];
rz(-1.7846196) q[1];
sx q[1];
rz(-2.7604483) q[1];
rz(-1.7140017) q[3];
sx q[3];
rz(-2.2369336) q[3];
sx q[3];
rz(-1.0216658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1537062) q[2];
sx q[2];
rz(-1.1517297) q[2];
sx q[2];
rz(0.64725867) q[2];
rz(-0.17793947) q[3];
sx q[3];
rz(-1.8837594) q[3];
sx q[3];
rz(1.2037163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7385638) q[0];
sx q[0];
rz(-0.084349923) q[0];
sx q[0];
rz(-0.52363288) q[0];
rz(1.9117484) q[1];
sx q[1];
rz(-0.74408999) q[1];
sx q[1];
rz(0.24807608) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6726674) q[0];
sx q[0];
rz(-1.4644196) q[0];
sx q[0];
rz(2.399978) q[0];
rz(-pi) q[1];
rz(-1.1311943) q[2];
sx q[2];
rz(-1.1710376) q[2];
sx q[2];
rz(-0.93610379) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9776002) q[1];
sx q[1];
rz(-1.5320141) q[1];
sx q[1];
rz(1.4250524) q[1];
rz(-2.4753597) q[3];
sx q[3];
rz(-1.1403475) q[3];
sx q[3];
rz(-2.7048049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.5292042) q[2];
sx q[2];
rz(-2.2233621) q[2];
sx q[2];
rz(0.53981346) q[2];
rz(-0.16048935) q[3];
sx q[3];
rz(-1.548998) q[3];
sx q[3];
rz(-0.81015712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6432583) q[0];
sx q[0];
rz(-2.46471) q[0];
sx q[0];
rz(-0.13370378) q[0];
rz(1.6224104) q[1];
sx q[1];
rz(-1.8459873) q[1];
sx q[1];
rz(-0.79230961) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2871029) q[0];
sx q[0];
rz(-0.40081319) q[0];
sx q[0];
rz(-2.0569747) q[0];
x q[1];
rz(2.4033264) q[2];
sx q[2];
rz(-0.19437899) q[2];
sx q[2];
rz(0.76972085) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6483602) q[1];
sx q[1];
rz(-0.91714749) q[1];
sx q[1];
rz(-2.3834989) q[1];
rz(-pi) q[2];
rz(0.00040690502) q[3];
sx q[3];
rz(-0.30012977) q[3];
sx q[3];
rz(0.89891922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2866659) q[2];
sx q[2];
rz(-2.4714405) q[2];
sx q[2];
rz(-0.83038846) q[2];
rz(-1.3487799) q[3];
sx q[3];
rz(-1.034779) q[3];
sx q[3];
rz(-2.5428298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8778359) q[0];
sx q[0];
rz(-0.97665518) q[0];
sx q[0];
rz(2.7184955) q[0];
rz(-1.9844203) q[1];
sx q[1];
rz(-1.6559699) q[1];
sx q[1];
rz(-0.62209904) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4870479) q[0];
sx q[0];
rz(-1.1830939) q[0];
sx q[0];
rz(0.88519208) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1272264) q[2];
sx q[2];
rz(-2.6831919) q[2];
sx q[2];
rz(1.0215525) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9378949) q[1];
sx q[1];
rz(-2.7337498) q[1];
sx q[1];
rz(1.6126584) q[1];
x q[2];
rz(-2.29633) q[3];
sx q[3];
rz(-1.5503395) q[3];
sx q[3];
rz(-1.6930228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.28079924) q[2];
sx q[2];
rz(-1.7223027) q[2];
sx q[2];
rz(0.61817509) q[2];
rz(1.6587967) q[3];
sx q[3];
rz(-1.7890472) q[3];
sx q[3];
rz(-1.6331204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16267714) q[0];
sx q[0];
rz(-1.8122883) q[0];
sx q[0];
rz(0.83258122) q[0];
rz(-2.816448) q[1];
sx q[1];
rz(-0.99601662) q[1];
sx q[1];
rz(3.0772193) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8125809) q[0];
sx q[0];
rz(-1.6906749) q[0];
sx q[0];
rz(2.8146312) q[0];
x q[1];
rz(-0.60466296) q[2];
sx q[2];
rz(-0.69398601) q[2];
sx q[2];
rz(2.1460905) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1015013) q[1];
sx q[1];
rz(-2.1143997) q[1];
sx q[1];
rz(0.60740791) q[1];
rz(-pi) q[2];
rz(-1.0945372) q[3];
sx q[3];
rz(-0.68449293) q[3];
sx q[3];
rz(1.5409926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.23318204) q[2];
sx q[2];
rz(-2.5711214) q[2];
sx q[2];
rz(1.2916267) q[2];
rz(2.6455247) q[3];
sx q[3];
rz(-0.78947624) q[3];
sx q[3];
rz(1.9991649) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18601501) q[0];
sx q[0];
rz(-2.8514974) q[0];
sx q[0];
rz(-2.7225323) q[0];
rz(2.8298607) q[1];
sx q[1];
rz(-1.5985951) q[1];
sx q[1];
rz(-0.14351235) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9547894) q[0];
sx q[0];
rz(-2.020103) q[0];
sx q[0];
rz(-1.9687998) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.2741998) q[2];
sx q[2];
rz(-1.3440709) q[2];
sx q[2];
rz(0.86967865) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4048481) q[1];
sx q[1];
rz(-0.84538922) q[1];
sx q[1];
rz(-2.0015697) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7147948) q[3];
sx q[3];
rz(-1.3132976) q[3];
sx q[3];
rz(-0.3443623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.36025563) q[2];
sx q[2];
rz(-0.90136734) q[2];
sx q[2];
rz(-2.0085013) q[2];
rz(1.1059149) q[3];
sx q[3];
rz(-1.4191041) q[3];
sx q[3];
rz(-2.6667986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.356242) q[0];
sx q[0];
rz(-0.79371912) q[0];
sx q[0];
rz(-1.445795) q[0];
rz(2.4244507) q[1];
sx q[1];
rz(-1.8627867) q[1];
sx q[1];
rz(0.76146567) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7324156) q[0];
sx q[0];
rz(-1.4237561) q[0];
sx q[0];
rz(0.63296403) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0324208) q[2];
sx q[2];
rz(-2.3512977) q[2];
sx q[2];
rz(-0.99144713) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.95132212) q[1];
sx q[1];
rz(-2.1441064) q[1];
sx q[1];
rz(0.058752937) q[1];
x q[2];
rz(-0.86628116) q[3];
sx q[3];
rz(-2.8968262) q[3];
sx q[3];
rz(0.94947986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.36934272) q[2];
sx q[2];
rz(-2.6476634) q[2];
sx q[2];
rz(0.93144766) q[2];
rz(-0.61819589) q[3];
sx q[3];
rz(-1.5074916) q[3];
sx q[3];
rz(0.96562323) q[3];
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
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6922927) q[0];
sx q[0];
rz(-1.3626784) q[0];
sx q[0];
rz(1.3344673) q[0];
rz(-2.6284699) q[1];
sx q[1];
rz(-2.158439) q[1];
sx q[1];
rz(1.2293053) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9216825) q[0];
sx q[0];
rz(-1.9905042) q[0];
sx q[0];
rz(-2.3786663) q[0];
rz(-pi) q[1];
rz(3.0487829) q[2];
sx q[2];
rz(-1.3121737) q[2];
sx q[2];
rz(1.5802204) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9521515) q[1];
sx q[1];
rz(-1.5501889) q[1];
sx q[1];
rz(2.9856696) q[1];
rz(0.66219285) q[3];
sx q[3];
rz(-1.5846888) q[3];
sx q[3];
rz(0.77869773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8339771) q[2];
sx q[2];
rz(-2.6145085) q[2];
sx q[2];
rz(0.21044593) q[2];
rz(-3.0757507) q[3];
sx q[3];
rz(-2.6688771) q[3];
sx q[3];
rz(-2.3426447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(0.54302067) q[0];
sx q[0];
rz(-0.6544756) q[0];
sx q[0];
rz(-1.2731113) q[0];
rz(-1.8674564) q[1];
sx q[1];
rz(-2.4576063) q[1];
sx q[1];
rz(-0.88409105) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40128517) q[0];
sx q[0];
rz(-1.570556) q[0];
sx q[0];
rz(-0.9539286) q[0];
x q[1];
rz(-2.7124518) q[2];
sx q[2];
rz(-2.57518) q[2];
sx q[2];
rz(-2.8265116) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6581304) q[1];
sx q[1];
rz(-2.6331055) q[1];
sx q[1];
rz(-2.1668651) q[1];
x q[2];
rz(0.69334778) q[3];
sx q[3];
rz(-2.7487488) q[3];
sx q[3];
rz(1.0885914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1511496) q[2];
sx q[2];
rz(-2.1572025) q[2];
sx q[2];
rz(1.3746064) q[2];
rz(-0.19045842) q[3];
sx q[3];
rz(-0.74646598) q[3];
sx q[3];
rz(1.9577352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9780438) q[0];
sx q[0];
rz(-1.4530285) q[0];
sx q[0];
rz(-1.2646041) q[0];
rz(-0.073607445) q[1];
sx q[1];
rz(-1.0818447) q[1];
sx q[1];
rz(0.61990613) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.093599565) q[0];
sx q[0];
rz(-1.7353829) q[0];
sx q[0];
rz(-0.82407804) q[0];
rz(1.0629547) q[2];
sx q[2];
rz(-1.0217102) q[2];
sx q[2];
rz(2.6428403) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3504847) q[1];
sx q[1];
rz(-1.3602166) q[1];
sx q[1];
rz(1.775022) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3315174) q[3];
sx q[3];
rz(-1.115723) q[3];
sx q[3];
rz(-0.51067715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8668883) q[2];
sx q[2];
rz(-0.70211774) q[2];
sx q[2];
rz(-2.9956024) q[2];
rz(-2.919096) q[3];
sx q[3];
rz(-0.22495088) q[3];
sx q[3];
rz(0.72993025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9324026) q[0];
sx q[0];
rz(-1.1840191) q[0];
sx q[0];
rz(-2.9970072) q[0];
rz(-1.2296386) q[1];
sx q[1];
rz(-1.3508136) q[1];
sx q[1];
rz(-1.2523686) q[1];
rz(-2.8242883) q[2];
sx q[2];
rz(-2.5226421) q[2];
sx q[2];
rz(0.68963827) q[2];
rz(1.7207809) q[3];
sx q[3];
rz(-0.21258988) q[3];
sx q[3];
rz(-1.7439738) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
