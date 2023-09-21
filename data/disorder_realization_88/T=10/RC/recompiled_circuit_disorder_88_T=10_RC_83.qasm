OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.8945583) q[0];
sx q[0];
rz(-1.2556827) q[0];
sx q[0];
rz(-0.32796252) q[0];
rz(-0.2344996) q[1];
sx q[1];
rz(-2.9405138) q[1];
sx q[1];
rz(-0.091436401) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2823921) q[0];
sx q[0];
rz(-0.99635591) q[0];
sx q[0];
rz(2.0496856) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.47317998) q[2];
sx q[2];
rz(-2.8521529) q[2];
sx q[2];
rz(-0.48130408) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4788471) q[1];
sx q[1];
rz(-1.4884243) q[1];
sx q[1];
rz(1.808951) q[1];
rz(-pi) q[2];
rz(-1.2535291) q[3];
sx q[3];
rz(-0.61875611) q[3];
sx q[3];
rz(1.3878824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4771007) q[2];
sx q[2];
rz(-2.1680919) q[2];
sx q[2];
rz(2.0155902) q[2];
rz(-2.8664355) q[3];
sx q[3];
rz(-2.5313009) q[3];
sx q[3];
rz(-0.90308213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1317516) q[0];
sx q[0];
rz(-2.4480902) q[0];
sx q[0];
rz(-2.4480208) q[0];
rz(1.0961078) q[1];
sx q[1];
rz(-0.98384905) q[1];
sx q[1];
rz(0.19031659) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41102558) q[0];
sx q[0];
rz(-1.1682296) q[0];
sx q[0];
rz(-0.21582027) q[0];
x q[1];
rz(1.9637738) q[2];
sx q[2];
rz(-1.0699546) q[2];
sx q[2];
rz(0.21265342) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9840934) q[1];
sx q[1];
rz(-1.5484973) q[1];
sx q[1];
rz(-0.84866546) q[1];
rz(-pi) q[2];
rz(2.6211561) q[3];
sx q[3];
rz(-2.3428829) q[3];
sx q[3];
rz(0.48898104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6341614) q[2];
sx q[2];
rz(-0.71175152) q[2];
sx q[2];
rz(1.2197536) q[2];
rz(-0.1427342) q[3];
sx q[3];
rz(-2.1295363) q[3];
sx q[3];
rz(-2.232961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74293566) q[0];
sx q[0];
rz(-2.0516899) q[0];
sx q[0];
rz(-0.8202585) q[0];
rz(0.29207686) q[1];
sx q[1];
rz(-2.0675817) q[1];
sx q[1];
rz(1.8935727) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4753715) q[0];
sx q[0];
rz(-0.9733805) q[0];
sx q[0];
rz(-1.682838) q[0];
x q[1];
rz(-2.9343611) q[2];
sx q[2];
rz(-1.5079632) q[2];
sx q[2];
rz(2.732423) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6229912) q[1];
sx q[1];
rz(-2.2366183) q[1];
sx q[1];
rz(1.4856505) q[1];
x q[2];
rz(1.1850584) q[3];
sx q[3];
rz(-2.4534561) q[3];
sx q[3];
rz(3.0582173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8212006) q[2];
sx q[2];
rz(-2.0911066) q[2];
sx q[2];
rz(-1.2134264) q[2];
rz(2.976867) q[3];
sx q[3];
rz(-2.2475593) q[3];
sx q[3];
rz(3.059982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7320025) q[0];
sx q[0];
rz(-1.859917) q[0];
sx q[0];
rz(0.18606342) q[0];
rz(-2.9371254) q[1];
sx q[1];
rz(-2.6696413) q[1];
sx q[1];
rz(1.2971372) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9545249) q[0];
sx q[0];
rz(-1.6141119) q[0];
sx q[0];
rz(-1.4403507) q[0];
rz(-pi) q[1];
x q[1];
rz(1.859971) q[2];
sx q[2];
rz(-1.7316069) q[2];
sx q[2];
rz(0.5493872) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0399196) q[1];
sx q[1];
rz(-2.1664201) q[1];
sx q[1];
rz(2.4458829) q[1];
x q[2];
rz(-1.7503731) q[3];
sx q[3];
rz(-1.5844126) q[3];
sx q[3];
rz(-2.7333958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.90594784) q[2];
sx q[2];
rz(-0.79259688) q[2];
sx q[2];
rz(0.91919351) q[2];
rz(2.8202608) q[3];
sx q[3];
rz(-2.0733757) q[3];
sx q[3];
rz(-1.2478158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9677251) q[0];
sx q[0];
rz(-0.50650948) q[0];
sx q[0];
rz(-2.2633973) q[0];
rz(1.325266) q[1];
sx q[1];
rz(-1.7293431) q[1];
sx q[1];
rz(1.7153046) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25709897) q[0];
sx q[0];
rz(-1.6528659) q[0];
sx q[0];
rz(-2.2280072) q[0];
x q[1];
rz(-2.2141453) q[2];
sx q[2];
rz(-2.7577835) q[2];
sx q[2];
rz(-0.62188934) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9702455) q[1];
sx q[1];
rz(-0.73927021) q[1];
sx q[1];
rz(1.7319748) q[1];
x q[2];
rz(2.812643) q[3];
sx q[3];
rz(-0.79862404) q[3];
sx q[3];
rz(0.22508276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2720126) q[2];
sx q[2];
rz(-11*pi/13) q[2];
sx q[2];
rz(-3.0997979) q[2];
rz(3.0801008) q[3];
sx q[3];
rz(-2.0525335) q[3];
sx q[3];
rz(2.7048236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7866661) q[0];
sx q[0];
rz(-3.0556398) q[0];
sx q[0];
rz(-2.1110995) q[0];
rz(-2.4018535) q[1];
sx q[1];
rz(-1.61295) q[1];
sx q[1];
rz(-2.5700263) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33486903) q[0];
sx q[0];
rz(-2.3248701) q[0];
sx q[0];
rz(-0.73210324) q[0];
rz(-1.7320485) q[2];
sx q[2];
rz(-0.67220062) q[2];
sx q[2];
rz(-1.6598998) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.9532277) q[1];
sx q[1];
rz(-1.1168915) q[1];
sx q[1];
rz(2.8737349) q[1];
rz(-pi) q[2];
rz(2.3267641) q[3];
sx q[3];
rz(-2.1216672) q[3];
sx q[3];
rz(-1.8967472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.968154) q[2];
sx q[2];
rz(-0.74192321) q[2];
sx q[2];
rz(-2.5773933) q[2];
rz(3.0155904) q[3];
sx q[3];
rz(-1.4583476) q[3];
sx q[3];
rz(2.8469767) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0512222) q[0];
sx q[0];
rz(-0.19989656) q[0];
sx q[0];
rz(-0.71227658) q[0];
rz(-2.6157216) q[1];
sx q[1];
rz(-2.7253175) q[1];
sx q[1];
rz(-0.66551048) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80488801) q[0];
sx q[0];
rz(-2.8993653) q[0];
sx q[0];
rz(1.5664943) q[0];
rz(-pi) q[1];
rz(-1.3940784) q[2];
sx q[2];
rz(-0.70509796) q[2];
sx q[2];
rz(1.9213898) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1674541) q[1];
sx q[1];
rz(-1.1145089) q[1];
sx q[1];
rz(2.8970701) q[1];
rz(-pi) q[2];
rz(-2.1036759) q[3];
sx q[3];
rz(-0.27927342) q[3];
sx q[3];
rz(-0.52475196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9843288) q[2];
sx q[2];
rz(-2.4815346) q[2];
sx q[2];
rz(-1.7187913) q[2];
rz(-0.11519365) q[3];
sx q[3];
rz(-2.6264103) q[3];
sx q[3];
rz(-2.9490024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.049906235) q[0];
sx q[0];
rz(-2.005907) q[0];
sx q[0];
rz(2.9507622) q[0];
rz(-2.514839) q[1];
sx q[1];
rz(-2.125506) q[1];
sx q[1];
rz(-2.802882) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7539464) q[0];
sx q[0];
rz(-0.37746143) q[0];
sx q[0];
rz(-1.1520555) q[0];
rz(0.45192265) q[2];
sx q[2];
rz(-1.3616614) q[2];
sx q[2];
rz(-2.506633) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0847816) q[1];
sx q[1];
rz(-1.4598795) q[1];
sx q[1];
rz(-1.0952428) q[1];
rz(-pi) q[2];
rz(1.2440153) q[3];
sx q[3];
rz(-1.8339001) q[3];
sx q[3];
rz(-0.14740482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.64615858) q[2];
sx q[2];
rz(-2.50714) q[2];
sx q[2];
rz(-0.70043606) q[2];
rz(2.2436079) q[3];
sx q[3];
rz(-1.286641) q[3];
sx q[3];
rz(0.17351304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53196466) q[0];
sx q[0];
rz(-0.48159972) q[0];
sx q[0];
rz(-0.65834808) q[0];
rz(0.61093962) q[1];
sx q[1];
rz(-1.2555723) q[1];
sx q[1];
rz(-0.13959612) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0848207) q[0];
sx q[0];
rz(-3.0773101) q[0];
sx q[0];
rz(-1.2540741) q[0];
rz(-pi) q[1];
x q[1];
rz(2.326968) q[2];
sx q[2];
rz(-1.6135718) q[2];
sx q[2];
rz(0.2864366) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.86233625) q[1];
sx q[1];
rz(-1.5039872) q[1];
sx q[1];
rz(-2.6194797) q[1];
rz(-pi) q[2];
rz(0.92630444) q[3];
sx q[3];
rz(-1.6127805) q[3];
sx q[3];
rz(-0.55461649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5247941) q[2];
sx q[2];
rz(-1.0107661) q[2];
sx q[2];
rz(-2.810478) q[2];
rz(0.75774276) q[3];
sx q[3];
rz(-2.7527633) q[3];
sx q[3];
rz(0.087879114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8452633) q[0];
sx q[0];
rz(-1.0110649) q[0];
sx q[0];
rz(2.9560126) q[0];
rz(1.0962076) q[1];
sx q[1];
rz(-2.9269693) q[1];
sx q[1];
rz(1.6569482) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.164924) q[0];
sx q[0];
rz(-0.52545588) q[0];
sx q[0];
rz(1.0111965) q[0];
rz(-2.2514258) q[2];
sx q[2];
rz(-2.7919127) q[2];
sx q[2];
rz(-2.4868951) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6757322) q[1];
sx q[1];
rz(-1.670174) q[1];
sx q[1];
rz(-1.194721) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5649892) q[3];
sx q[3];
rz(-0.70492893) q[3];
sx q[3];
rz(0.73202902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.93402702) q[2];
sx q[2];
rz(-2.2521314) q[2];
sx q[2];
rz(-0.55220848) q[2];
rz(-2.3637555) q[3];
sx q[3];
rz(-2.2838897) q[3];
sx q[3];
rz(-2.623693) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26375297) q[0];
sx q[0];
rz(-1.629566) q[0];
sx q[0];
rz(-1.4019479) q[0];
rz(0.18763018) q[1];
sx q[1];
rz(-1.7788806) q[1];
sx q[1];
rz(2.3685041) q[1];
rz(-0.27992237) q[2];
sx q[2];
rz(-1.7114077) q[2];
sx q[2];
rz(-0.68243295) q[2];
rz(-2.256176) q[3];
sx q[3];
rz(-2.1366742) q[3];
sx q[3];
rz(-2.4921806) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];