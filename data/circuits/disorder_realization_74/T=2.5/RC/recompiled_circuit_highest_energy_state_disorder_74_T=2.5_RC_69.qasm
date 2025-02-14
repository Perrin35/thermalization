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
rz(0.49630961) q[1];
sx q[1];
rz(-0.30975431) q[1];
sx q[1];
rz(0.53909477) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2236299) q[0];
sx q[0];
rz(-1.5962036) q[0];
sx q[0];
rz(-0.44141234) q[0];
rz(-pi) q[1];
rz(-2.7819958) q[2];
sx q[2];
rz(-1.7126377) q[2];
sx q[2];
rz(-1.8388909) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8546363) q[1];
sx q[1];
rz(-2.7883734) q[1];
sx q[1];
rz(-1.1562721) q[1];
x q[2];
rz(-2.21978) q[3];
sx q[3];
rz(-1.6361437) q[3];
sx q[3];
rz(-1.1999038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3794136) q[2];
sx q[2];
rz(-2.2654686) q[2];
sx q[2];
rz(1.1890821) q[2];
rz(-1.1842229) q[3];
sx q[3];
rz(-0.90457478) q[3];
sx q[3];
rz(-1.5949465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.9963843) q[0];
sx q[0];
rz(-1.5897911) q[0];
sx q[0];
rz(1.8183964) q[0];
rz(-2.6595751) q[1];
sx q[1];
rz(-0.91164416) q[1];
sx q[1];
rz(2.1673896) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8708987) q[0];
sx q[0];
rz(-2.2091731) q[0];
sx q[0];
rz(-2.9186072) q[0];
rz(3.0423099) q[2];
sx q[2];
rz(-1.6834604) q[2];
sx q[2];
rz(2.7275865) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5017101) q[1];
sx q[1];
rz(-0.5941662) q[1];
sx q[1];
rz(-0.42826786) q[1];
rz(0.5663381) q[3];
sx q[3];
rz(-1.7741331) q[3];
sx q[3];
rz(-1.7254064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.0049858967) q[2];
sx q[2];
rz(-0.58711457) q[2];
sx q[2];
rz(2.3919487) q[2];
rz(0.43198112) q[3];
sx q[3];
rz(-1.1155198) q[3];
sx q[3];
rz(-1.3031134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62464803) q[0];
sx q[0];
rz(-0.2722781) q[0];
sx q[0];
rz(-0.6063478) q[0];
rz(1.8079405) q[1];
sx q[1];
rz(-1.9275815) q[1];
sx q[1];
rz(2.8820754) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5657046) q[0];
sx q[0];
rz(-1.4153) q[0];
sx q[0];
rz(1.3208273) q[0];
x q[1];
rz(-2.9664842) q[2];
sx q[2];
rz(-1.7705743) q[2];
sx q[2];
rz(-1.4005043) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.77713457) q[1];
sx q[1];
rz(-1.196047) q[1];
sx q[1];
rz(1.115429) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7068038) q[3];
sx q[3];
rz(-2.090095) q[3];
sx q[3];
rz(0.10552191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2567265) q[2];
sx q[2];
rz(-1.3639516) q[2];
sx q[2];
rz(-1.5474896) q[2];
rz(0.78440845) q[3];
sx q[3];
rz(-1.6268077) q[3];
sx q[3];
rz(1.5659531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6467658) q[0];
sx q[0];
rz(-2.0443199) q[0];
sx q[0];
rz(0.25949091) q[0];
rz(1.2207813) q[1];
sx q[1];
rz(-0.44993284) q[1];
sx q[1];
rz(1.6993914) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9913015) q[0];
sx q[0];
rz(-1.1027005) q[0];
sx q[0];
rz(2.1670114) q[0];
x q[1];
rz(-3.0109826) q[2];
sx q[2];
rz(-2.161986) q[2];
sx q[2];
rz(-1.3863877) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.759023) q[1];
sx q[1];
rz(-1.8174531) q[1];
sx q[1];
rz(-0.36212977) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.070017858) q[3];
sx q[3];
rz(-2.2517831) q[3];
sx q[3];
rz(0.79047608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1059025) q[2];
sx q[2];
rz(-1.8041958) q[2];
sx q[2];
rz(2.7316459) q[2];
rz(-1.9299054) q[3];
sx q[3];
rz(-2.4544921) q[3];
sx q[3];
rz(1.80779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7164417) q[0];
sx q[0];
rz(-0.71617675) q[0];
sx q[0];
rz(2.8675365) q[0];
rz(2.7033499) q[1];
sx q[1];
rz(-1.848105) q[1];
sx q[1];
rz(-0.73572198) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.05984027) q[0];
sx q[0];
rz(-0.81938374) q[0];
sx q[0];
rz(-1.6900586) q[0];
rz(-1.7154791) q[2];
sx q[2];
rz(-1.4229454) q[2];
sx q[2];
rz(1.3593909) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3408518) q[1];
sx q[1];
rz(-1.2781118) q[1];
sx q[1];
rz(-0.26009788) q[1];
rz(-2.6295119) q[3];
sx q[3];
rz(-2.1554073) q[3];
sx q[3];
rz(-1.1325815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2030486) q[2];
sx q[2];
rz(-1.6837348) q[2];
sx q[2];
rz(2.5277444) q[2];
rz(-0.40890536) q[3];
sx q[3];
rz(-0.96858612) q[3];
sx q[3];
rz(-0.74234211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3530389) q[0];
sx q[0];
rz(-2.5034294) q[0];
sx q[0];
rz(1.2370538) q[0];
rz(1.6963814) q[1];
sx q[1];
rz(-0.45672363) q[1];
sx q[1];
rz(2.9663185) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4458081) q[0];
sx q[0];
rz(-0.17477594) q[0];
sx q[0];
rz(1.4265027) q[0];
x q[1];
rz(-2.1138328) q[2];
sx q[2];
rz(-0.98688417) q[2];
sx q[2];
rz(1.8915382) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1239667) q[1];
sx q[1];
rz(-2.0793036) q[1];
sx q[1];
rz(2.070049) q[1];
rz(-pi) q[2];
x q[2];
rz(0.82476576) q[3];
sx q[3];
rz(-2.3769393) q[3];
sx q[3];
rz(-0.95461707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.96482977) q[2];
sx q[2];
rz(-1.333933) q[2];
sx q[2];
rz(-0.97243398) q[2];
rz(1.4771627) q[3];
sx q[3];
rz(-1.9921314) q[3];
sx q[3];
rz(0.76178637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2870188) q[0];
sx q[0];
rz(-0.022495689) q[0];
sx q[0];
rz(0.35312411) q[0];
rz(-2.1137721) q[1];
sx q[1];
rz(-1.2007583) q[1];
sx q[1];
rz(1.0110528) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42732692) q[0];
sx q[0];
rz(-1.2135226) q[0];
sx q[0];
rz(1.5861804) q[0];
x q[1];
rz(0.89247668) q[2];
sx q[2];
rz(-2.3111642) q[2];
sx q[2];
rz(2.3811946) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.331703) q[1];
sx q[1];
rz(-1.2223489) q[1];
sx q[1];
rz(2.9747559) q[1];
rz(-pi) q[2];
rz(0.8259186) q[3];
sx q[3];
rz(-2.7248235) q[3];
sx q[3];
rz(0.59405223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8280243) q[2];
sx q[2];
rz(-2.821533) q[2];
sx q[2];
rz(1.774452) q[2];
rz(1.8732871) q[3];
sx q[3];
rz(-1.6627848) q[3];
sx q[3];
rz(0.84793276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1139514) q[0];
sx q[0];
rz(-2.5340762) q[0];
sx q[0];
rz(-1.129958) q[0];
rz(0.21895151) q[1];
sx q[1];
rz(-1.7349225) q[1];
sx q[1];
rz(0.91046441) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.509678) q[0];
sx q[0];
rz(-0.70075894) q[0];
sx q[0];
rz(0.52978446) q[0];
x q[1];
rz(2.6633419) q[2];
sx q[2];
rz(-1.646981) q[2];
sx q[2];
rz(-2.8752799) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1840802) q[1];
sx q[1];
rz(-1.0614479) q[1];
sx q[1];
rz(-0.30708509) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.72515709) q[3];
sx q[3];
rz(-1.2598383) q[3];
sx q[3];
rz(-2.6127315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3160481) q[2];
sx q[2];
rz(-2.3393708) q[2];
sx q[2];
rz(2.3160589) q[2];
rz(0.46142203) q[3];
sx q[3];
rz(-1.2666707) q[3];
sx q[3];
rz(0.44574827) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
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
rz(-0.83183944) q[0];
rz(0.72744751) q[1];
sx q[1];
rz(-1.4201545) q[1];
sx q[1];
rz(0.096253455) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.03503) q[0];
sx q[0];
rz(-2.0756531) q[0];
sx q[0];
rz(-0.52953203) q[0];
rz(-0.25419323) q[2];
sx q[2];
rz(-1.5113748) q[2];
sx q[2];
rz(-2.3583902) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0301275) q[1];
sx q[1];
rz(-0.22997738) q[1];
sx q[1];
rz(1.1632009) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7670125) q[3];
sx q[3];
rz(-1.726578) q[3];
sx q[3];
rz(2.2957612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7353797) q[2];
sx q[2];
rz(-1.4344183) q[2];
sx q[2];
rz(1.5926788) q[2];
rz(2.5088572) q[3];
sx q[3];
rz(-1.9097208) q[3];
sx q[3];
rz(-0.24937853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-11/(9*pi)) q[0];
sx q[0];
rz(-1.0405552) q[0];
sx q[0];
rz(2.7885875) q[0];
rz(2.8189335) q[1];
sx q[1];
rz(-2.8162075) q[1];
sx q[1];
rz(-0.37193146) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3561365) q[0];
sx q[0];
rz(-1.8710941) q[0];
sx q[0];
rz(2.9343283) q[0];
rz(2.6549321) q[2];
sx q[2];
rz(-2.7722205) q[2];
sx q[2];
rz(-0.83115679) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5185753) q[1];
sx q[1];
rz(-1.2378344) q[1];
sx q[1];
rz(1.740231) q[1];
rz(-pi) q[2];
rz(-0.84623611) q[3];
sx q[3];
rz(-1.6000984) q[3];
sx q[3];
rz(-1.8054363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9670664) q[2];
sx q[2];
rz(-1.9517978) q[2];
sx q[2];
rz(1.2971499) q[2];
rz(-0.8477115) q[3];
sx q[3];
rz(-1.1573557) q[3];
sx q[3];
rz(-2.1367836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-0.90947718) q[0];
sx q[0];
rz(-1.7625325) q[0];
sx q[0];
rz(0.43176227) q[0];
rz(-1.3879981) q[1];
sx q[1];
rz(-1.2139865) q[1];
sx q[1];
rz(-0.075275631) q[1];
rz(-0.83609381) q[2];
sx q[2];
rz(-0.89955624) q[2];
sx q[2];
rz(1.6513165) q[2];
rz(-1.9907436) q[3];
sx q[3];
rz(-2.3227878) q[3];
sx q[3];
rz(-0.25521758) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
