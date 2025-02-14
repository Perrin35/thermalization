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
rz(0.81646252) q[0];
sx q[0];
rz(-3.0397968) q[0];
sx q[0];
rz(0.53959227) q[0];
rz(3.6379023) q[1];
sx q[1];
rz(3.451347) q[1];
sx q[1];
rz(5.7440905) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2236299) q[0];
sx q[0];
rz(-1.5453891) q[0];
sx q[0];
rz(-0.44141234) q[0];
x q[1];
rz(1.7221872) q[2];
sx q[2];
rz(-1.2149723) q[2];
sx q[2];
rz(-0.21499888) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8546363) q[1];
sx q[1];
rz(-0.35321924) q[1];
sx q[1];
rz(-1.9853206) q[1];
rz(-pi) q[2];
rz(-2.21978) q[3];
sx q[3];
rz(-1.505449) q[3];
sx q[3];
rz(-1.9416888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3794136) q[2];
sx q[2];
rz(-0.87612408) q[2];
sx q[2];
rz(-1.1890821) q[2];
rz(-1.9573697) q[3];
sx q[3];
rz(-0.90457478) q[3];
sx q[3];
rz(-1.5466461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14520833) q[0];
sx q[0];
rz(-1.5518016) q[0];
sx q[0];
rz(-1.3231963) q[0];
rz(0.48201758) q[1];
sx q[1];
rz(-0.91164416) q[1];
sx q[1];
rz(-0.97420305) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4344113) q[0];
sx q[0];
rz(-1.392258) q[0];
sx q[0];
rz(-0.92037998) q[0];
rz(-pi) q[1];
rz(3.0423099) q[2];
sx q[2];
rz(-1.6834604) q[2];
sx q[2];
rz(2.7275865) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5017101) q[1];
sx q[1];
rz(-0.5941662) q[1];
sx q[1];
rz(2.7133248) q[1];
rz(-pi) q[2];
rz(1.810435) q[3];
sx q[3];
rz(-2.1240892) q[3];
sx q[3];
rz(-3.1146793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.1366068) q[2];
sx q[2];
rz(-2.5544781) q[2];
sx q[2];
rz(-0.74964398) q[2];
rz(2.7096115) q[3];
sx q[3];
rz(-2.0260729) q[3];
sx q[3];
rz(1.8384793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62464803) q[0];
sx q[0];
rz(-0.2722781) q[0];
sx q[0];
rz(0.6063478) q[0];
rz(1.8079405) q[1];
sx q[1];
rz(-1.9275815) q[1];
sx q[1];
rz(-0.25951728) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6919428) q[0];
sx q[0];
rz(-2.8480673) q[0];
sx q[0];
rz(-2.13563) q[0];
rz(-0.86028966) q[2];
sx q[2];
rz(-0.26488556) q[2];
sx q[2];
rz(2.4691894) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.77713457) q[1];
sx q[1];
rz(-1.9455457) q[1];
sx q[1];
rz(1.115429) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0084148) q[3];
sx q[3];
rz(-1.1964238) q[3];
sx q[3];
rz(-1.9028456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2567265) q[2];
sx q[2];
rz(-1.3639516) q[2];
sx q[2];
rz(-1.5941031) q[2];
rz(2.3571842) q[3];
sx q[3];
rz(-1.514785) q[3];
sx q[3];
rz(-1.5756395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6467658) q[0];
sx q[0];
rz(-2.0443199) q[0];
sx q[0];
rz(0.25949091) q[0];
rz(-1.2207813) q[1];
sx q[1];
rz(-0.44993284) q[1];
sx q[1];
rz(1.4422013) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15029112) q[0];
sx q[0];
rz(-1.1027005) q[0];
sx q[0];
rz(-2.1670114) q[0];
x q[1];
rz(-2.1659508) q[2];
sx q[2];
rz(-1.6791428) q[2];
sx q[2];
rz(0.25749046) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3842786) q[1];
sx q[1];
rz(-0.43506778) q[1];
sx q[1];
rz(-0.61788322) q[1];
rz(-2.2529834) q[3];
sx q[3];
rz(-1.6251792) q[3];
sx q[3];
rz(0.82444421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1059025) q[2];
sx q[2];
rz(-1.3373969) q[2];
sx q[2];
rz(0.40994677) q[2];
rz(1.9299054) q[3];
sx q[3];
rz(-0.68710059) q[3];
sx q[3];
rz(1.80779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7164417) q[0];
sx q[0];
rz(-0.71617675) q[0];
sx q[0];
rz(2.8675365) q[0];
rz(-2.7033499) q[1];
sx q[1];
rz(-1.848105) q[1];
sx q[1];
rz(0.73572198) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0817524) q[0];
sx q[0];
rz(-0.81938374) q[0];
sx q[0];
rz(1.6900586) q[0];
x q[1];
rz(-2.9922036) q[2];
sx q[2];
rz(-1.7138897) q[2];
sx q[2];
rz(-2.9516475) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.80074089) q[1];
sx q[1];
rz(-1.8634808) q[1];
sx q[1];
rz(2.8814948) q[1];
x q[2];
rz(0.93345668) q[3];
sx q[3];
rz(-2.3846755) q[3];
sx q[3];
rz(2.8035823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.938544) q[2];
sx q[2];
rz(-1.6837348) q[2];
sx q[2];
rz(2.5277444) q[2];
rz(-0.40890536) q[3];
sx q[3];
rz(-2.1730065) q[3];
sx q[3];
rz(-2.3992505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3530389) q[0];
sx q[0];
rz(-2.5034294) q[0];
sx q[0];
rz(-1.9045389) q[0];
rz(1.4452112) q[1];
sx q[1];
rz(-2.684869) q[1];
sx q[1];
rz(-0.17527418) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4458081) q[0];
sx q[0];
rz(-2.9668167) q[0];
sx q[0];
rz(-1.4265027) q[0];
x q[1];
rz(2.1138328) q[2];
sx q[2];
rz(-2.1547085) q[2];
sx q[2];
rz(1.8915382) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8126909) q[1];
sx q[1];
rz(-2.0022503) q[1];
sx q[1];
rz(-0.56568362) q[1];
x q[2];
rz(-0.95703362) q[3];
sx q[3];
rz(-2.0599457) q[3];
sx q[3];
rz(-0.027994284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.96482977) q[2];
sx q[2];
rz(-1.333933) q[2];
sx q[2];
rz(0.97243398) q[2];
rz(1.6644299) q[3];
sx q[3];
rz(-1.1494613) q[3];
sx q[3];
rz(-2.3798063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85457388) q[0];
sx q[0];
rz(-0.022495689) q[0];
sx q[0];
rz(-2.7884685) q[0];
rz(1.0278206) q[1];
sx q[1];
rz(-1.2007583) q[1];
sx q[1];
rz(1.0110528) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38336223) q[0];
sx q[0];
rz(-0.35759059) q[0];
sx q[0];
rz(3.1004058) q[0];
rz(-pi) q[1];
rz(-0.6017466) q[2];
sx q[2];
rz(-2.1831552) q[2];
sx q[2];
rz(-1.5073206) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8450635) q[1];
sx q[1];
rz(-1.7275212) q[1];
sx q[1];
rz(1.2178376) q[1];
x q[2];
rz(1.8854463) q[3];
sx q[3];
rz(-1.2928179) q[3];
sx q[3];
rz(-1.6772456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8280243) q[2];
sx q[2];
rz(-0.32005969) q[2];
sx q[2];
rz(1.3671406) q[2];
rz(1.8732871) q[3];
sx q[3];
rz(-1.6627848) q[3];
sx q[3];
rz(-2.2936599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0276412) q[0];
sx q[0];
rz(-2.5340762) q[0];
sx q[0];
rz(-1.129958) q[0];
rz(2.9226411) q[1];
sx q[1];
rz(-1.7349225) q[1];
sx q[1];
rz(2.2311282) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2856871) q[0];
sx q[0];
rz(-2.160851) q[0];
sx q[0];
rz(1.1678215) q[0];
rz(-pi) q[1];
rz(-1.4850281) q[2];
sx q[2];
rz(-1.0940486) q[2];
sx q[2];
rz(-1.7976744) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1840802) q[1];
sx q[1];
rz(-1.0614479) q[1];
sx q[1];
rz(0.30708509) q[1];
rz(-0.72515709) q[3];
sx q[3];
rz(-1.8817543) q[3];
sx q[3];
rz(2.6127315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.8255446) q[2];
sx q[2];
rz(-2.3393708) q[2];
sx q[2];
rz(0.82553378) q[2];
rz(2.6801706) q[3];
sx q[3];
rz(-1.2666707) q[3];
sx q[3];
rz(2.6958444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6398741) q[0];
sx q[0];
rz(-1.3383144) q[0];
sx q[0];
rz(2.3097532) q[0];
rz(0.72744751) q[1];
sx q[1];
rz(-1.7214382) q[1];
sx q[1];
rz(-0.096253455) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1883586) q[0];
sx q[0];
rz(-2.0287477) q[0];
sx q[0];
rz(1.0012549) q[0];
x q[1];
rz(0.25419323) q[2];
sx q[2];
rz(-1.5113748) q[2];
sx q[2];
rz(2.3583902) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.93864307) q[1];
sx q[1];
rz(-1.4803107) q[1];
sx q[1];
rz(1.7825148) q[1];
rz(-pi) q[2];
rz(1.3745802) q[3];
sx q[3];
rz(-1.726578) q[3];
sx q[3];
rz(2.2957612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7353797) q[2];
sx q[2];
rz(-1.4344183) q[2];
sx q[2];
rz(-1.5489138) q[2];
rz(-0.63273543) q[3];
sx q[3];
rz(-1.2318719) q[3];
sx q[3];
rz(-2.8922141) q[3];
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
rz(2.7525472) q[0];
sx q[0];
rz(-1.0405552) q[0];
sx q[0];
rz(2.7885875) q[0];
rz(0.32265916) q[1];
sx q[1];
rz(-0.32538515) q[1];
sx q[1];
rz(2.7696612) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8648099) q[0];
sx q[0];
rz(-1.3729334) q[0];
sx q[0];
rz(1.2643306) q[0];
rz(0.48666059) q[2];
sx q[2];
rz(-0.36937215) q[2];
sx q[2];
rz(2.3104359) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5185753) q[1];
sx q[1];
rz(-1.2378344) q[1];
sx q[1];
rz(-1.4013616) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6149893) q[3];
sx q[3];
rz(-2.4165476) q[3];
sx q[3];
rz(-2.8738662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9670664) q[2];
sx q[2];
rz(-1.9517978) q[2];
sx q[2];
rz(-1.8444427) q[2];
rz(-2.2938812) q[3];
sx q[3];
rz(-1.984237) q[3];
sx q[3];
rz(1.004809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90947718) q[0];
sx q[0];
rz(-1.7625325) q[0];
sx q[0];
rz(0.43176227) q[0];
rz(1.7535946) q[1];
sx q[1];
rz(-1.2139865) q[1];
sx q[1];
rz(-0.075275631) q[1];
rz(2.3222011) q[2];
sx q[2];
rz(-1.0181622) q[2];
sx q[2];
rz(-0.43140585) q[2];
rz(-2.344178) q[3];
sx q[3];
rz(-1.8731464) q[3];
sx q[3];
rz(-2.1220589) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
