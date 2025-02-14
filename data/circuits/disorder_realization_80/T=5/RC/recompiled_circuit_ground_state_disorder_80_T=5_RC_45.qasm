OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.6516946) q[0];
sx q[0];
rz(-0.88180056) q[0];
sx q[0];
rz(-0.14468004) q[0];
rz(0.41809234) q[1];
sx q[1];
rz(-0.10377181) q[1];
sx q[1];
rz(1.511908) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9448175) q[0];
sx q[0];
rz(-0.51048764) q[0];
sx q[0];
rz(1.8148242) q[0];
rz(-pi) q[1];
rz(-1.51654) q[2];
sx q[2];
rz(-2.3035604) q[2];
sx q[2];
rz(-2.4861479) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.92951951) q[1];
sx q[1];
rz(-1.1843425) q[1];
sx q[1];
rz(-2.7824336) q[1];
rz(2.6295794) q[3];
sx q[3];
rz(-1.9422934) q[3];
sx q[3];
rz(-0.12646261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8251553) q[2];
sx q[2];
rz(-1.8708159) q[2];
sx q[2];
rz(-0.6553418) q[2];
rz(1.3842899) q[3];
sx q[3];
rz(-0.96450788) q[3];
sx q[3];
rz(2.3206319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4493988) q[0];
sx q[0];
rz(-1.7658424) q[0];
sx q[0];
rz(-2.6334515) q[0];
rz(-1.3805768) q[1];
sx q[1];
rz(-2.5602129) q[1];
sx q[1];
rz(1.2095721) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0469477) q[0];
sx q[0];
rz(-2.1535465) q[0];
sx q[0];
rz(2.9215982) q[0];
x q[1];
rz(-0.31866535) q[2];
sx q[2];
rz(-0.056431596) q[2];
sx q[2];
rz(2.3669764) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3928814) q[1];
sx q[1];
rz(-1.6986753) q[1];
sx q[1];
rz(-0.46462469) q[1];
rz(-pi) q[2];
rz(2.1749635) q[3];
sx q[3];
rz(-1.894891) q[3];
sx q[3];
rz(1.5965727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9380583) q[2];
sx q[2];
rz(-1.7455696) q[2];
sx q[2];
rz(-2.5246942) q[2];
rz(1.6054224) q[3];
sx q[3];
rz(-1.0441531) q[3];
sx q[3];
rz(-0.49083403) q[3];
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
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9817552) q[0];
sx q[0];
rz(-1.6438537) q[0];
sx q[0];
rz(-2.4165261) q[0];
rz(1.2695351) q[1];
sx q[1];
rz(-0.67765647) q[1];
sx q[1];
rz(-2.1824172) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3676783) q[0];
sx q[0];
rz(-1.1836021) q[0];
sx q[0];
rz(-2.0810803) q[0];
rz(2.3940635) q[2];
sx q[2];
rz(-0.5420712) q[2];
sx q[2];
rz(0.67491787) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.37466418) q[1];
sx q[1];
rz(-2.9247724) q[1];
sx q[1];
rz(-2.1310852) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7531583) q[3];
sx q[3];
rz(-1.7373457) q[3];
sx q[3];
rz(2.9349851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8150669) q[2];
sx q[2];
rz(-1.8748137) q[2];
sx q[2];
rz(-0.25700021) q[2];
rz(1.0810931) q[3];
sx q[3];
rz(-0.5210146) q[3];
sx q[3];
rz(1.5628373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0839888) q[0];
sx q[0];
rz(-2.8631518) q[0];
sx q[0];
rz(1.9637015) q[0];
rz(0.15011694) q[1];
sx q[1];
rz(-2.6094486) q[1];
sx q[1];
rz(1.6252801) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1193084) q[0];
sx q[0];
rz(-1.2890576) q[0];
sx q[0];
rz(-1.9153829) q[0];
x q[1];
rz(-0.84211911) q[2];
sx q[2];
rz(-1.2884022) q[2];
sx q[2];
rz(-0.51674622) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.1399556) q[1];
sx q[1];
rz(-1.6703342) q[1];
sx q[1];
rz(2.8147354) q[1];
x q[2];
rz(2.4978906) q[3];
sx q[3];
rz(-1.4182404) q[3];
sx q[3];
rz(-1.2920472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.84245044) q[2];
sx q[2];
rz(-2.5717042) q[2];
sx q[2];
rz(2.2184856) q[2];
rz(-2.182377) q[3];
sx q[3];
rz(-1.0126637) q[3];
sx q[3];
rz(2.925351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8206896) q[0];
sx q[0];
rz(-2.237759) q[0];
sx q[0];
rz(-0.88291105) q[0];
rz(1.8461022) q[1];
sx q[1];
rz(-2.3675282) q[1];
sx q[1];
rz(-0.49497089) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8710324) q[0];
sx q[0];
rz(-2.705276) q[0];
sx q[0];
rz(-2.3113234) q[0];
rz(-2.1033435) q[2];
sx q[2];
rz(-0.94013273) q[2];
sx q[2];
rz(-1.3884461) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.043083) q[1];
sx q[1];
rz(-1.6585697) q[1];
sx q[1];
rz(0.66941525) q[1];
rz(-pi) q[2];
x q[2];
rz(0.99397387) q[3];
sx q[3];
rz(-0.60189542) q[3];
sx q[3];
rz(-2.8771597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3317269) q[2];
sx q[2];
rz(-0.11955424) q[2];
sx q[2];
rz(1.2681819) q[2];
rz(3.0590893) q[3];
sx q[3];
rz(-1.5039597) q[3];
sx q[3];
rz(-0.65703195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6562011) q[0];
sx q[0];
rz(-2.9347561) q[0];
sx q[0];
rz(-1.6400826) q[0];
rz(-1.062475) q[1];
sx q[1];
rz(-1.1524009) q[1];
sx q[1];
rz(1.9662439) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15099111) q[0];
sx q[0];
rz(-2.7620533) q[0];
sx q[0];
rz(-2.5113456) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7926867) q[2];
sx q[2];
rz(-2.4173173) q[2];
sx q[2];
rz(-1.5021715) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3835241) q[1];
sx q[1];
rz(-1.822377) q[1];
sx q[1];
rz(1.4497767) q[1];
x q[2];
rz(-2.8296521) q[3];
sx q[3];
rz(-0.7850724) q[3];
sx q[3];
rz(-1.5842445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8106653) q[2];
sx q[2];
rz(-1.284548) q[2];
sx q[2];
rz(-0.98768273) q[2];
rz(-2.2641613) q[3];
sx q[3];
rz(-2.8743447) q[3];
sx q[3];
rz(-2.4115244) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1838609) q[0];
sx q[0];
rz(-1.6022302) q[0];
sx q[0];
rz(0.13885942) q[0];
rz(-1.214437) q[1];
sx q[1];
rz(-1.5077533) q[1];
sx q[1];
rz(0.81659281) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79838419) q[0];
sx q[0];
rz(-1.3460165) q[0];
sx q[0];
rz(-0.019545743) q[0];
rz(-2.8726878) q[2];
sx q[2];
rz(-2.3550526) q[2];
sx q[2];
rz(1.5603017) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.066799435) q[1];
sx q[1];
rz(-1.6082676) q[1];
sx q[1];
rz(2.2897359) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4754627) q[3];
sx q[3];
rz(-0.96530789) q[3];
sx q[3];
rz(0.44119409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.82178086) q[2];
sx q[2];
rz(-2.767441) q[2];
sx q[2];
rz(0.4846586) q[2];
rz(-2.7759077) q[3];
sx q[3];
rz(-1.7771143) q[3];
sx q[3];
rz(1.1588089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0115688) q[0];
sx q[0];
rz(-2.8304709) q[0];
sx q[0];
rz(-3.044361) q[0];
rz(-0.10487996) q[1];
sx q[1];
rz(-0.77969867) q[1];
sx q[1];
rz(1.3146776) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86634462) q[0];
sx q[0];
rz(-1.5684248) q[0];
sx q[0];
rz(1.5644685) q[0];
rz(1.8123367) q[2];
sx q[2];
rz(-1.353423) q[2];
sx q[2];
rz(-1.5940983) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.842031) q[1];
sx q[1];
rz(-0.52134575) q[1];
sx q[1];
rz(3.035665) q[1];
rz(-3.0544326) q[3];
sx q[3];
rz(-1.0356552) q[3];
sx q[3];
rz(-2.5453107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9647727) q[2];
sx q[2];
rz(-1.7440045) q[2];
sx q[2];
rz(-2.9883265) q[2];
rz(-1.2166474) q[3];
sx q[3];
rz(-0.21806589) q[3];
sx q[3];
rz(2.5254068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.1185054) q[0];
sx q[0];
rz(-3.1361134) q[0];
sx q[0];
rz(1.5047005) q[0];
rz(-2.3432689) q[1];
sx q[1];
rz(-1.6183805) q[1];
sx q[1];
rz(-0.33531478) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42683168) q[0];
sx q[0];
rz(-2.8726467) q[0];
sx q[0];
rz(1.8105276) q[0];
x q[1];
rz(1.5952571) q[2];
sx q[2];
rz(-0.82055295) q[2];
sx q[2];
rz(-0.94708196) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0287231) q[1];
sx q[1];
rz(-1.2838893) q[1];
sx q[1];
rz(2.7863281) q[1];
rz(-pi) q[2];
rz(-0.9934523) q[3];
sx q[3];
rz(-1.7029188) q[3];
sx q[3];
rz(-0.41228774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.12751427) q[2];
sx q[2];
rz(-1.9731584) q[2];
sx q[2];
rz(0.43759313) q[2];
rz(-2.0994999) q[3];
sx q[3];
rz(-1.4053248) q[3];
sx q[3];
rz(2.5991345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9893148) q[0];
sx q[0];
rz(-2.2570026) q[0];
sx q[0];
rz(0.44961318) q[0];
rz(0.97688976) q[1];
sx q[1];
rz(-1.5039109) q[1];
sx q[1];
rz(-0.50122112) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64912632) q[0];
sx q[0];
rz(-1.7607795) q[0];
sx q[0];
rz(1.282592) q[0];
x q[1];
rz(1.1468723) q[2];
sx q[2];
rz(-1.3739283) q[2];
sx q[2];
rz(-2.9483861) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.75647351) q[1];
sx q[1];
rz(-1.3344527) q[1];
sx q[1];
rz(2.0467351) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.95285033) q[3];
sx q[3];
rz(-2.5898159) q[3];
sx q[3];
rz(-2.1976041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0235128) q[2];
sx q[2];
rz(-2.2492275) q[2];
sx q[2];
rz(-0.21772131) q[2];
rz(-1.9067541) q[3];
sx q[3];
rz(-0.66399884) q[3];
sx q[3];
rz(0.92065221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6476743) q[0];
sx q[0];
rz(-1.7128581) q[0];
sx q[0];
rz(0.38059522) q[0];
rz(-2.4299798) q[1];
sx q[1];
rz(-1.8743534) q[1];
sx q[1];
rz(1.4030917) q[1];
rz(0.099945036) q[2];
sx q[2];
rz(-2.6700085) q[2];
sx q[2];
rz(2.811583) q[2];
rz(-1.3344516) q[3];
sx q[3];
rz(-1.9665039) q[3];
sx q[3];
rz(2.891091) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
