OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6310298) q[0];
sx q[0];
rz(-2.6369542) q[0];
sx q[0];
rz(-2.7161427) q[0];
rz(0.57234859) q[1];
sx q[1];
rz(4.2387716) q[1];
sx q[1];
rz(7.4833202) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10353032) q[0];
sx q[0];
rz(-1.5804884) q[0];
sx q[0];
rz(-0.25083019) q[0];
rz(2.8699371) q[2];
sx q[2];
rz(-1.2363529) q[2];
sx q[2];
rz(-1.7890499) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7537743) q[1];
sx q[1];
rz(-1.942831) q[1];
sx q[1];
rz(-1.340999) q[1];
x q[2];
rz(-0.17961924) q[3];
sx q[3];
rz(-2.4625375) q[3];
sx q[3];
rz(-1.8906901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9878865) q[2];
sx q[2];
rz(-1.989863) q[2];
sx q[2];
rz(-2.494334) q[2];
rz(0.17793947) q[3];
sx q[3];
rz(-1.2578332) q[3];
sx q[3];
rz(-1.9378763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7385638) q[0];
sx q[0];
rz(-0.084349923) q[0];
sx q[0];
rz(0.52363288) q[0];
rz(1.2298443) q[1];
sx q[1];
rz(-2.3975027) q[1];
sx q[1];
rz(-2.8935166) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46892525) q[0];
sx q[0];
rz(-1.6771731) q[0];
sx q[0];
rz(-0.74161462) q[0];
rz(-2.0103983) q[2];
sx q[2];
rz(-1.9705551) q[2];
sx q[2];
rz(-0.93610379) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9776002) q[1];
sx q[1];
rz(-1.5320141) q[1];
sx q[1];
rz(1.4250524) q[1];
rz(-pi) q[2];
rz(2.0994151) q[3];
sx q[3];
rz(-0.97441593) q[3];
sx q[3];
rz(-1.6906052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6123885) q[2];
sx q[2];
rz(-2.2233621) q[2];
sx q[2];
rz(2.6017792) q[2];
rz(0.16048935) q[3];
sx q[3];
rz(-1.548998) q[3];
sx q[3];
rz(0.81015712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4983343) q[0];
sx q[0];
rz(-2.46471) q[0];
sx q[0];
rz(0.13370378) q[0];
rz(-1.6224104) q[1];
sx q[1];
rz(-1.8459873) q[1];
sx q[1];
rz(-2.349283) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8544898) q[0];
sx q[0];
rz(-0.40081319) q[0];
sx q[0];
rz(-2.0569747) q[0];
rz(-pi) q[1];
rz(-2.4033264) q[2];
sx q[2];
rz(-0.19437899) q[2];
sx q[2];
rz(2.3718718) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.50666675) q[1];
sx q[1];
rz(-2.185195) q[1];
sx q[1];
rz(-2.3022815) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5706704) q[3];
sx q[3];
rz(-1.8709261) q[3];
sx q[3];
rz(-2.2430994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2866659) q[2];
sx q[2];
rz(-2.4714405) q[2];
sx q[2];
rz(2.3112042) q[2];
rz(-1.3487799) q[3];
sx q[3];
rz(-2.1068137) q[3];
sx q[3];
rz(2.5428298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26375672) q[0];
sx q[0];
rz(-0.97665518) q[0];
sx q[0];
rz(0.42309716) q[0];
rz(1.1571723) q[1];
sx q[1];
rz(-1.4856228) q[1];
sx q[1];
rz(-2.5194936) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34919993) q[0];
sx q[0];
rz(-2.3696941) q[0];
sx q[0];
rz(-0.99794879) q[0];
x q[1];
rz(1.1515205) q[2];
sx q[2];
rz(-1.3797233) q[2];
sx q[2];
rz(2.9951468) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8129261) q[1];
sx q[1];
rz(-1.587396) q[1];
sx q[1];
rz(-1.9783201) q[1];
rz(-pi) q[2];
rz(0.84526269) q[3];
sx q[3];
rz(-1.5912531) q[3];
sx q[3];
rz(-1.4485698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.28079924) q[2];
sx q[2];
rz(-1.7223027) q[2];
sx q[2];
rz(-2.5234176) q[2];
rz(-1.6587967) q[3];
sx q[3];
rz(-1.7890472) q[3];
sx q[3];
rz(1.6331204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16267714) q[0];
sx q[0];
rz(-1.8122883) q[0];
sx q[0];
rz(-2.3090114) q[0];
rz(-2.816448) q[1];
sx q[1];
rz(-0.99601662) q[1];
sx q[1];
rz(3.0772193) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8125809) q[0];
sx q[0];
rz(-1.4509177) q[0];
sx q[0];
rz(0.32696149) q[0];
rz(-pi) q[1];
rz(-0.60466296) q[2];
sx q[2];
rz(-0.69398601) q[2];
sx q[2];
rz(-0.99550216) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0325614) q[1];
sx q[1];
rz(-2.3500725) q[1];
sx q[1];
rz(2.3276106) q[1];
rz(-pi) q[2];
rz(1.0945372) q[3];
sx q[3];
rz(-0.68449293) q[3];
sx q[3];
rz(1.6006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.23318204) q[2];
sx q[2];
rz(-2.5711214) q[2];
sx q[2];
rz(1.849966) q[2];
rz(0.49606797) q[3];
sx q[3];
rz(-2.3521164) q[3];
sx q[3];
rz(-1.1424278) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18601501) q[0];
sx q[0];
rz(-2.8514974) q[0];
sx q[0];
rz(-2.7225323) q[0];
rz(0.31173197) q[1];
sx q[1];
rz(-1.5985951) q[1];
sx q[1];
rz(0.14351235) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3449609) q[0];
sx q[0];
rz(-1.9274184) q[0];
sx q[0];
rz(-0.48194569) q[0];
x q[1];
rz(-2.8673929) q[2];
sx q[2];
rz(-1.7975217) q[2];
sx q[2];
rz(0.86967865) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.12998768) q[1];
sx q[1];
rz(-1.2530348) q[1];
sx q[1];
rz(-0.77316534) q[1];
rz(-2.6426598) q[3];
sx q[3];
rz(-0.29424516) q[3];
sx q[3];
rz(0.1732451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.781337) q[2];
sx q[2];
rz(-0.90136734) q[2];
sx q[2];
rz(2.0085013) q[2];
rz(-2.0356778) q[3];
sx q[3];
rz(-1.4191041) q[3];
sx q[3];
rz(0.47479409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7853506) q[0];
sx q[0];
rz(-0.79371912) q[0];
sx q[0];
rz(-1.6957977) q[0];
rz(-2.4244507) q[1];
sx q[1];
rz(-1.278806) q[1];
sx q[1];
rz(-2.380127) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7324156) q[0];
sx q[0];
rz(-1.7178365) q[0];
sx q[0];
rz(0.63296403) q[0];
rz(0.85651969) q[2];
sx q[2];
rz(-1.9437143) q[2];
sx q[2];
rz(-0.97717092) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1902705) q[1];
sx q[1];
rz(-2.1441064) q[1];
sx q[1];
rz(0.058752937) q[1];
rz(-2.2753115) q[3];
sx q[3];
rz(-2.8968262) q[3];
sx q[3];
rz(-0.94947986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.36934272) q[2];
sx q[2];
rz(-0.49392924) q[2];
sx q[2];
rz(0.93144766) q[2];
rz(0.61819589) q[3];
sx q[3];
rz(-1.5074916) q[3];
sx q[3];
rz(2.1759694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6922927) q[0];
sx q[0];
rz(-1.3626784) q[0];
sx q[0];
rz(-1.8071254) q[0];
rz(-0.5131228) q[1];
sx q[1];
rz(-0.98315364) q[1];
sx q[1];
rz(-1.9122874) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4192159) q[0];
sx q[0];
rz(-0.88787365) q[0];
sx q[0];
rz(-1.0177259) q[0];
rz(-pi) q[1];
rz(1.3111054) q[2];
sx q[2];
rz(-1.4810815) q[2];
sx q[2];
rz(-0.033223991) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1894412) q[1];
sx q[1];
rz(-1.5914038) q[1];
sx q[1];
rz(-0.15592305) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.66219285) q[3];
sx q[3];
rz(-1.5846888) q[3];
sx q[3];
rz(2.3628949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3076155) q[2];
sx q[2];
rz(-2.6145085) q[2];
sx q[2];
rz(2.9311467) q[2];
rz(3.0757507) q[3];
sx q[3];
rz(-0.47271553) q[3];
sx q[3];
rz(0.79894799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.598572) q[0];
sx q[0];
rz(-2.4871171) q[0];
sx q[0];
rz(1.2731113) q[0];
rz(-1.2741362) q[1];
sx q[1];
rz(-0.68398634) q[1];
sx q[1];
rz(-0.88409105) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1696816) q[0];
sx q[0];
rz(-2.187664) q[0];
sx q[0];
rz(3.141298) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.52427788) q[2];
sx q[2];
rz(-1.7959716) q[2];
sx q[2];
rz(0.88722992) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6581304) q[1];
sx q[1];
rz(-0.50848714) q[1];
sx q[1];
rz(-2.1668651) q[1];
x q[2];
rz(0.69334778) q[3];
sx q[3];
rz(-0.39284387) q[3];
sx q[3];
rz(-1.0885914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.990443) q[2];
sx q[2];
rz(-2.1572025) q[2];
sx q[2];
rz(-1.7669862) q[2];
rz(-2.9511342) q[3];
sx q[3];
rz(-2.3951267) q[3];
sx q[3];
rz(1.9577352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16354887) q[0];
sx q[0];
rz(-1.4530285) q[0];
sx q[0];
rz(1.8769886) q[0];
rz(0.073607445) q[1];
sx q[1];
rz(-2.059748) q[1];
sx q[1];
rz(0.61990613) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0479931) q[0];
sx q[0];
rz(-1.7353829) q[0];
sx q[0];
rz(-0.82407804) q[0];
rz(-0.61087278) q[2];
sx q[2];
rz(-1.1430642) q[2];
sx q[2];
rz(-1.3547225) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3504847) q[1];
sx q[1];
rz(-1.7813761) q[1];
sx q[1];
rz(1.775022) q[1];
x q[2];
rz(-2.6750071) q[3];
sx q[3];
rz(-1.3562725) q[3];
sx q[3];
rz(1.1669351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8668883) q[2];
sx q[2];
rz(-0.70211774) q[2];
sx q[2];
rz(-0.14599027) q[2];
rz(0.22249666) q[3];
sx q[3];
rz(-0.22495088) q[3];
sx q[3];
rz(-2.4116624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20919007) q[0];
sx q[0];
rz(-1.9575735) q[0];
sx q[0];
rz(0.14458543) q[0];
rz(1.9119541) q[1];
sx q[1];
rz(-1.3508136) q[1];
sx q[1];
rz(-1.2523686) q[1];
rz(0.59496224) q[2];
sx q[2];
rz(-1.7528201) q[2];
sx q[2];
rz(2.5217944) q[2];
rz(0.03224198) q[3];
sx q[3];
rz(-1.7809636) q[3];
sx q[3];
rz(-1.8973593) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
