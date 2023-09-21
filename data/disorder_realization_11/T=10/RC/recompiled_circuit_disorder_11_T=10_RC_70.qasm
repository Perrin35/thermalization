OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.55460632) q[0];
sx q[0];
rz(-0.89590961) q[0];
sx q[0];
rz(1.7370261) q[0];
rz(0.70146927) q[1];
sx q[1];
rz(-0.62086064) q[1];
sx q[1];
rz(1.4863185) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0044999997) q[0];
sx q[0];
rz(-0.21919964) q[0];
sx q[0];
rz(2.5704434) q[0];
rz(-2.9831224) q[2];
sx q[2];
rz(-2.7140693) q[2];
sx q[2];
rz(2.1536364) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.15581407) q[1];
sx q[1];
rz(-1.474838) q[1];
sx q[1];
rz(-0.39246651) q[1];
rz(-1.5641883) q[3];
sx q[3];
rz(-1.4918543) q[3];
sx q[3];
rz(3.1066432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0554589) q[2];
sx q[2];
rz(-1.3749342) q[2];
sx q[2];
rz(1.3936183) q[2];
rz(0.80111516) q[3];
sx q[3];
rz(-1.5603742) q[3];
sx q[3];
rz(-2.0786409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3515117) q[0];
sx q[0];
rz(-1.693049) q[0];
sx q[0];
rz(-3.063607) q[0];
rz(0.69411913) q[1];
sx q[1];
rz(-1.1342987) q[1];
sx q[1];
rz(2.7819113) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65950899) q[0];
sx q[0];
rz(-0.88920278) q[0];
sx q[0];
rz(-0.73007749) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4215135) q[2];
sx q[2];
rz(-1.4125925) q[2];
sx q[2];
rz(0.0062696487) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.82275326) q[1];
sx q[1];
rz(-2.0944632) q[1];
sx q[1];
rz(0.19284064) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.85675591) q[3];
sx q[3];
rz(-1.5407011) q[3];
sx q[3];
rz(1.4853958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2461207) q[2];
sx q[2];
rz(-1.3161696) q[2];
sx q[2];
rz(0.94334156) q[2];
rz(-1.881276) q[3];
sx q[3];
rz(-2.3337119) q[3];
sx q[3];
rz(2.4356306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27186069) q[0];
sx q[0];
rz(-1.6534709) q[0];
sx q[0];
rz(0.46762064) q[0];
rz(-0.46472654) q[1];
sx q[1];
rz(-0.48370353) q[1];
sx q[1];
rz(2.8288249) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5178459) q[0];
sx q[0];
rz(-0.22652921) q[0];
sx q[0];
rz(1.9342599) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8094693) q[2];
sx q[2];
rz(-2.2687074) q[2];
sx q[2];
rz(0.78813625) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.73156563) q[1];
sx q[1];
rz(-1.7419852) q[1];
sx q[1];
rz(-2.9982135) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7160277) q[3];
sx q[3];
rz(-1.8947453) q[3];
sx q[3];
rz(-2.5026929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.23190817) q[2];
sx q[2];
rz(-1.568305) q[2];
sx q[2];
rz(2.8364733) q[2];
rz(-1.8049847) q[3];
sx q[3];
rz(-0.87953416) q[3];
sx q[3];
rz(0.46521503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6453648) q[0];
sx q[0];
rz(-2.7858758) q[0];
sx q[0];
rz(-0.028045068) q[0];
rz(1.6758945) q[1];
sx q[1];
rz(-2.5391255) q[1];
sx q[1];
rz(2.4688597) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69986491) q[0];
sx q[0];
rz(-1.2727951) q[0];
sx q[0];
rz(2.5900048) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5920934) q[2];
sx q[2];
rz(-1.3504343) q[2];
sx q[2];
rz(-2.9133177) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.87186253) q[1];
sx q[1];
rz(-1.6849663) q[1];
sx q[1];
rz(-1.8037379) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1898515) q[3];
sx q[3];
rz(-1.674106) q[3];
sx q[3];
rz(2.2593474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4611886) q[2];
sx q[2];
rz(-2.350312) q[2];
sx q[2];
rz(2.380774) q[2];
rz(-0.86756724) q[3];
sx q[3];
rz(-1.8224199) q[3];
sx q[3];
rz(-2.3755465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5714394) q[0];
sx q[0];
rz(-1.1076936) q[0];
sx q[0];
rz(-1.0636348) q[0];
rz(-1.4798374) q[1];
sx q[1];
rz(-0.96264833) q[1];
sx q[1];
rz(3.1033049) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83513658) q[0];
sx q[0];
rz(-1.4493754) q[0];
sx q[0];
rz(1.2760713) q[0];
x q[1];
rz(-2.4942336) q[2];
sx q[2];
rz(-0.72509662) q[2];
sx q[2];
rz(2.4908096) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2513411) q[1];
sx q[1];
rz(-2.9311935) q[1];
sx q[1];
rz(0.9743147) q[1];
rz(2.8713545) q[3];
sx q[3];
rz(-1.1186557) q[3];
sx q[3];
rz(-0.42735344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.15497196) q[2];
sx q[2];
rz(-1.9382696) q[2];
sx q[2];
rz(-1.0920452) q[2];
rz(1.5345705) q[3];
sx q[3];
rz(-2.1708596) q[3];
sx q[3];
rz(0.78305125) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18096481) q[0];
sx q[0];
rz(-1.6406849) q[0];
sx q[0];
rz(-1.3076179) q[0];
rz(0.062782137) q[1];
sx q[1];
rz(-1.9654013) q[1];
sx q[1];
rz(-2.9506156) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4673429) q[0];
sx q[0];
rz(-1.4313587) q[0];
sx q[0];
rz(-1.3583202) q[0];
rz(-pi) q[1];
rz(0.88626185) q[2];
sx q[2];
rz(-1.2949416) q[2];
sx q[2];
rz(-0.72826284) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.2497017) q[1];
sx q[1];
rz(-0.54348031) q[1];
sx q[1];
rz(1.6103423) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1359547) q[3];
sx q[3];
rz(-1.8668381) q[3];
sx q[3];
rz(-2.5067096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6993616) q[2];
sx q[2];
rz(-1.8403534) q[2];
sx q[2];
rz(0.43506518) q[2];
rz(2.2502031) q[3];
sx q[3];
rz(-1.9458709) q[3];
sx q[3];
rz(-1.5589327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1973535) q[0];
sx q[0];
rz(-2.9159912) q[0];
sx q[0];
rz(-2.7213668) q[0];
rz(0.48121437) q[1];
sx q[1];
rz(-2.2599506) q[1];
sx q[1];
rz(-1.0302672) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97354613) q[0];
sx q[0];
rz(-1.9424115) q[0];
sx q[0];
rz(-0.51460534) q[0];
rz(-1.3517802) q[2];
sx q[2];
rz(-1.4957301) q[2];
sx q[2];
rz(-0.53596562) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.87228862) q[1];
sx q[1];
rz(-1.8708806) q[1];
sx q[1];
rz(-1.3809204) q[1];
rz(-pi) q[2];
rz(-1.7610735) q[3];
sx q[3];
rz(-1.7576302) q[3];
sx q[3];
rz(-2.2639757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.39945012) q[2];
sx q[2];
rz(-2.7059677) q[2];
sx q[2];
rz(0.31526652) q[2];
rz(2.1049843) q[3];
sx q[3];
rz(-1.2549812) q[3];
sx q[3];
rz(2.172327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2675562) q[0];
sx q[0];
rz(-0.75796217) q[0];
sx q[0];
rz(-1.1918921) q[0];
rz(-2.2299956) q[1];
sx q[1];
rz(-2.4688768) q[1];
sx q[1];
rz(-1.0620767) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9174791) q[0];
sx q[0];
rz(-2.5878721) q[0];
sx q[0];
rz(-1.8138767) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9757134) q[2];
sx q[2];
rz(-1.9654044) q[2];
sx q[2];
rz(-1.6369866) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.83795611) q[1];
sx q[1];
rz(-1.0567697) q[1];
sx q[1];
rz(2.282826) q[1];
x q[2];
rz(1.6154556) q[3];
sx q[3];
rz(-2.4756845) q[3];
sx q[3];
rz(-0.33223104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.186211) q[2];
sx q[2];
rz(-2.8583128) q[2];
sx q[2];
rz(-3.0905159) q[2];
rz(-0.71930277) q[3];
sx q[3];
rz(-1.5161113) q[3];
sx q[3];
rz(-2.0843845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99701571) q[0];
sx q[0];
rz(-0.41467312) q[0];
sx q[0];
rz(2.4966519) q[0];
rz(-2.0058517) q[1];
sx q[1];
rz(-0.93808162) q[1];
sx q[1];
rz(-0.84916806) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97201559) q[0];
sx q[0];
rz(-1.9398085) q[0];
sx q[0];
rz(-0.24032648) q[0];
rz(-pi) q[1];
rz(-1.4184985) q[2];
sx q[2];
rz(-1.2010152) q[2];
sx q[2];
rz(2.5068482) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4918684) q[1];
sx q[1];
rz(-2.8102311) q[1];
sx q[1];
rz(-2.3359873) q[1];
x q[2];
rz(-2.455515) q[3];
sx q[3];
rz(-0.25732532) q[3];
sx q[3];
rz(-2.0285332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7944305) q[2];
sx q[2];
rz(-2.590245) q[2];
sx q[2];
rz(0.17865044) q[2];
rz(-0.84154877) q[3];
sx q[3];
rz(-2.0984564) q[3];
sx q[3];
rz(-2.6508022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0586044) q[0];
sx q[0];
rz(-1.937979) q[0];
sx q[0];
rz(-0.9151181) q[0];
rz(-1.8944342) q[1];
sx q[1];
rz(-1.9688537) q[1];
sx q[1];
rz(-2.9291005) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1157163) q[0];
sx q[0];
rz(-1.3256782) q[0];
sx q[0];
rz(-2.2239767) q[0];
rz(-1.882952) q[2];
sx q[2];
rz(-1.6098795) q[2];
sx q[2];
rz(2.8651819) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.64582981) q[1];
sx q[1];
rz(-2.1730575) q[1];
sx q[1];
rz(-2.5219003) q[1];
rz(-1.8260164) q[3];
sx q[3];
rz(-0.85382429) q[3];
sx q[3];
rz(-2.9361801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.96182573) q[2];
sx q[2];
rz(-0.13040725) q[2];
sx q[2];
rz(1.4546222) q[2];
rz(-1.9466594) q[3];
sx q[3];
rz(-1.8353728) q[3];
sx q[3];
rz(-0.52836829) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4029978) q[0];
sx q[0];
rz(-1.4017372) q[0];
sx q[0];
rz(-1.5177939) q[0];
rz(3.0659499) q[1];
sx q[1];
rz(-1.6041926) q[1];
sx q[1];
rz(1.4354482) q[1];
rz(3.1235789) q[2];
sx q[2];
rz(-1.5009038) q[2];
sx q[2];
rz(0.76138087) q[2];
rz(0.71804071) q[3];
sx q[3];
rz(-2.4951039) q[3];
sx q[3];
rz(-1.390355) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];