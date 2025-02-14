OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4898981) q[0];
sx q[0];
rz(-2.2597921) q[0];
sx q[0];
rz(0.14468004) q[0];
rz(0.41809234) q[1];
sx q[1];
rz(-0.10377181) q[1];
sx q[1];
rz(1.511908) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83990012) q[0];
sx q[0];
rz(-1.6891251) q[0];
sx q[0];
rz(2.0685643) q[0];
rz(-1.51654) q[2];
sx q[2];
rz(-0.83803229) q[2];
sx q[2];
rz(2.4861479) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.92951951) q[1];
sx q[1];
rz(-1.1843425) q[1];
sx q[1];
rz(2.7824336) q[1];
x q[2];
rz(-1.9910664) q[3];
sx q[3];
rz(-2.0448677) q[3];
sx q[3];
rz(1.898511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.31643733) q[2];
sx q[2];
rz(-1.2707767) q[2];
sx q[2];
rz(2.4862508) q[2];
rz(-1.3842899) q[3];
sx q[3];
rz(-2.1770848) q[3];
sx q[3];
rz(2.3206319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4493988) q[0];
sx q[0];
rz(-1.7658424) q[0];
sx q[0];
rz(-0.50814116) q[0];
rz(-1.3805768) q[1];
sx q[1];
rz(-2.5602129) q[1];
sx q[1];
rz(1.2095721) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0469477) q[0];
sx q[0];
rz(-0.98804615) q[0];
sx q[0];
rz(-0.21999448) q[0];
x q[1];
rz(2.8229273) q[2];
sx q[2];
rz(-0.056431596) q[2];
sx q[2];
rz(2.3669764) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0703549) q[1];
sx q[1];
rz(-2.6609328) q[1];
sx q[1];
rz(2.8621469) q[1];
x q[2];
rz(-0.38755667) q[3];
sx q[3];
rz(-1.0021375) q[3];
sx q[3];
rz(2.89944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9380583) q[2];
sx q[2];
rz(-1.396023) q[2];
sx q[2];
rz(2.5246942) q[2];
rz(1.5361702) q[3];
sx q[3];
rz(-2.0974396) q[3];
sx q[3];
rz(2.6507586) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9817552) q[0];
sx q[0];
rz(-1.497739) q[0];
sx q[0];
rz(0.72506654) q[0];
rz(-1.8720576) q[1];
sx q[1];
rz(-0.67765647) q[1];
sx q[1];
rz(-2.1824172) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1364258) q[0];
sx q[0];
rz(-1.1014897) q[0];
sx q[0];
rz(0.43715663) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.959402) q[2];
sx q[2];
rz(-1.1827785) q[2];
sx q[2];
rz(-1.4997375) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.19634464) q[1];
sx q[1];
rz(-1.7540534) q[1];
sx q[1];
rz(0.11653479) q[1];
x q[2];
rz(0.16930468) q[3];
sx q[3];
rz(-1.3909855) q[3];
sx q[3];
rz(-1.7468417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8150669) q[2];
sx q[2];
rz(-1.2667789) q[2];
sx q[2];
rz(2.8845924) q[2];
rz(-2.0604996) q[3];
sx q[3];
rz(-0.5210146) q[3];
sx q[3];
rz(-1.5787554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0576039) q[0];
sx q[0];
rz(-0.27844089) q[0];
sx q[0];
rz(-1.1778911) q[0];
rz(-2.9914757) q[1];
sx q[1];
rz(-2.6094486) q[1];
sx q[1];
rz(1.6252801) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0310136) q[0];
sx q[0];
rz(-2.7000801) q[0];
sx q[0];
rz(2.2792321) q[0];
x q[1];
rz(2.2994735) q[2];
sx q[2];
rz(-1.8531905) q[2];
sx q[2];
rz(0.51674622) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8575688) q[1];
sx q[1];
rz(-2.800436) q[1];
sx q[1];
rz(-2.8400303) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4978906) q[3];
sx q[3];
rz(-1.4182404) q[3];
sx q[3];
rz(1.8495454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.84245044) q[2];
sx q[2];
rz(-2.5717042) q[2];
sx q[2];
rz(-2.2184856) q[2];
rz(-2.182377) q[3];
sx q[3];
rz(-1.0126637) q[3];
sx q[3];
rz(-0.21624163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8206896) q[0];
sx q[0];
rz(-2.237759) q[0];
sx q[0];
rz(2.2586816) q[0];
rz(-1.8461022) q[1];
sx q[1];
rz(-2.3675282) q[1];
sx q[1];
rz(-2.6466218) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51909101) q[0];
sx q[0];
rz(-1.8880196) q[0];
sx q[0];
rz(2.836801) q[0];
rz(1.0382492) q[2];
sx q[2];
rz(-2.2014599) q[2];
sx q[2];
rz(1.3884461) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.40303142) q[1];
sx q[1];
rz(-0.90442362) q[1];
sx q[1];
rz(-1.4590461) q[1];
rz(-pi) q[2];
rz(-0.35843973) q[3];
sx q[3];
rz(-2.0652986) q[3];
sx q[3];
rz(-2.2089437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8098658) q[2];
sx q[2];
rz(-0.11955424) q[2];
sx q[2];
rz(1.2681819) q[2];
rz(-0.082503334) q[3];
sx q[3];
rz(-1.5039597) q[3];
sx q[3];
rz(2.4845607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6562011) q[0];
sx q[0];
rz(-0.20683658) q[0];
sx q[0];
rz(-1.6400826) q[0];
rz(2.0791176) q[1];
sx q[1];
rz(-1.9891918) q[1];
sx q[1];
rz(1.1753488) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1262681) q[0];
sx q[0];
rz(-1.3506753) q[0];
sx q[0];
rz(0.31173978) q[0];
rz(-1.2771036) q[2];
sx q[2];
rz(-0.8986462) q[2];
sx q[2];
rz(-1.1872684) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.7580686) q[1];
sx q[1];
rz(-1.3192156) q[1];
sx q[1];
rz(1.4497767) q[1];
x q[2];
rz(-0.76035108) q[3];
sx q[3];
rz(-1.3521128) q[3];
sx q[3];
rz(-2.9038871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8106653) q[2];
sx q[2];
rz(-1.8570447) q[2];
sx q[2];
rz(-2.1539099) q[2];
rz(2.2641613) q[3];
sx q[3];
rz(-0.26724795) q[3];
sx q[3];
rz(0.73006829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9577318) q[0];
sx q[0];
rz(-1.6022302) q[0];
sx q[0];
rz(0.13885942) q[0];
rz(-1.9271556) q[1];
sx q[1];
rz(-1.6338394) q[1];
sx q[1];
rz(0.81659281) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3432085) q[0];
sx q[0];
rz(-1.7955762) q[0];
sx q[0];
rz(0.019545743) q[0];
rz(-pi) q[1];
rz(1.8310407) q[2];
sx q[2];
rz(-0.81965551) q[2];
sx q[2];
rz(-1.2092839) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6703709) q[1];
sx q[1];
rz(-2.2891217) q[1];
sx q[1];
rz(3.0918151) q[1];
x q[2];
rz(2.4754627) q[3];
sx q[3];
rz(-2.1762848) q[3];
sx q[3];
rz(2.7003986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3198118) q[2];
sx q[2];
rz(-2.767441) q[2];
sx q[2];
rz(2.6569341) q[2];
rz(2.7759077) q[3];
sx q[3];
rz(-1.7771143) q[3];
sx q[3];
rz(-1.1588089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1300238) q[0];
sx q[0];
rz(-0.31112177) q[0];
sx q[0];
rz(3.044361) q[0];
rz(-0.10487996) q[1];
sx q[1];
rz(-2.361894) q[1];
sx q[1];
rz(1.8269151) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0630232) q[0];
sx q[0];
rz(-0.0067575909) q[0];
sx q[0];
rz(1.2122173) q[0];
rz(0.22366053) q[2];
sx q[2];
rz(-1.8065435) q[2];
sx q[2];
rz(3.0652114) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7200249) q[1];
sx q[1];
rz(-2.0889258) q[1];
sx q[1];
rz(-1.5101449) q[1];
rz(1.7165887) q[3];
sx q[3];
rz(-2.6000826) q[3];
sx q[3];
rz(-2.7150115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9647727) q[2];
sx q[2];
rz(-1.7440045) q[2];
sx q[2];
rz(-2.9883265) q[2];
rz(1.9249453) q[3];
sx q[3];
rz(-2.9235268) q[3];
sx q[3];
rz(0.6161859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1185054) q[0];
sx q[0];
rz(-0.0054792976) q[0];
sx q[0];
rz(-1.5047005) q[0];
rz(0.79832375) q[1];
sx q[1];
rz(-1.6183805) q[1];
sx q[1];
rz(-0.33531478) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3753877) q[0];
sx q[0];
rz(-1.5076625) q[0];
sx q[0];
rz(1.8324018) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1153572) q[2];
sx q[2];
rz(-0.75056428) q[2];
sx q[2];
rz(-2.158643) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7882725) q[1];
sx q[1];
rz(-1.9109374) q[1];
sx q[1];
rz(1.8756833) q[1];
rz(-pi) q[2];
rz(1.8096381) q[3];
sx q[3];
rz(-0.59058978) q[3];
sx q[3];
rz(2.1826377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0140784) q[2];
sx q[2];
rz(-1.1684343) q[2];
sx q[2];
rz(2.7039995) q[2];
rz(-1.0420927) q[3];
sx q[3];
rz(-1.4053248) q[3];
sx q[3];
rz(-2.5991345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9893148) q[0];
sx q[0];
rz(-0.88459009) q[0];
sx q[0];
rz(0.44961318) q[0];
rz(2.1647029) q[1];
sx q[1];
rz(-1.6376817) q[1];
sx q[1];
rz(2.6403715) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86574449) q[0];
sx q[0];
rz(-1.8536708) q[0];
sx q[0];
rz(-2.9436443) q[0];
rz(-pi) q[1];
rz(1.1193163) q[2];
sx q[2];
rz(-0.46487936) q[2];
sx q[2];
rz(-1.7865739) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.75647351) q[1];
sx q[1];
rz(-1.8071399) q[1];
sx q[1];
rz(-1.0948576) q[1];
rz(-0.34256713) q[3];
sx q[3];
rz(-1.1293355) q[3];
sx q[3];
rz(1.502069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1180798) q[2];
sx q[2];
rz(-0.89236516) q[2];
sx q[2];
rz(-2.9238713) q[2];
rz(-1.2348385) q[3];
sx q[3];
rz(-2.4775938) q[3];
sx q[3];
rz(0.92065221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
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
rz(2.4299798) q[1];
sx q[1];
rz(-1.2672392) q[1];
sx q[1];
rz(-1.738501) q[1];
rz(-3.0416476) q[2];
sx q[2];
rz(-2.6700085) q[2];
sx q[2];
rz(2.811583) q[2];
rz(2.7357581) q[3];
sx q[3];
rz(-1.7885699) q[3];
sx q[3];
rz(-1.7287265) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
