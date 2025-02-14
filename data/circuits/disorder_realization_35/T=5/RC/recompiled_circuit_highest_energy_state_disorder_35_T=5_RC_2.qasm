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
rz(1.9995243) q[1];
sx q[1];
rz(-2.6231782) q[1];
sx q[1];
rz(1.4919182) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60580743) q[0];
sx q[0];
rz(-1.5174179) q[0];
sx q[0];
rz(-3.1230631) q[0];
rz(-pi) q[1];
rz(-2.1989735) q[2];
sx q[2];
rz(-0.97731263) q[2];
sx q[2];
rz(1.817199) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8644997) q[1];
sx q[1];
rz(-2.1756209) q[1];
sx q[1];
rz(-0.34326174) q[1];
rz(-pi) q[2];
rz(-3.1077216) q[3];
sx q[3];
rz(-1.813214) q[3];
sx q[3];
rz(-2.9592379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7923183) q[2];
sx q[2];
rz(-0.26733843) q[2];
sx q[2];
rz(0.054923687) q[2];
rz(1.4548291) q[3];
sx q[3];
rz(-0.39595404) q[3];
sx q[3];
rz(1.5168064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9216264) q[0];
sx q[0];
rz(-1.8333789) q[0];
sx q[0];
rz(2.8935905) q[0];
rz(-1.8964881) q[1];
sx q[1];
rz(-0.73687941) q[1];
sx q[1];
rz(-1.9128333) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22511521) q[0];
sx q[0];
rz(-0.18813293) q[0];
sx q[0];
rz(-2.9636901) q[0];
rz(-0.80557786) q[2];
sx q[2];
rz(-2.4586939) q[2];
sx q[2];
rz(2.3685092) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.23763021) q[1];
sx q[1];
rz(-1.7824689) q[1];
sx q[1];
rz(-2.4304076) q[1];
x q[2];
rz(-0.42028479) q[3];
sx q[3];
rz(-1.5160504) q[3];
sx q[3];
rz(-1.960567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9241141) q[2];
sx q[2];
rz(-0.83196297) q[2];
sx q[2];
rz(0.58383101) q[2];
rz(-2.7122688) q[3];
sx q[3];
rz(-1.2238294) q[3];
sx q[3];
rz(2.1302285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9871224) q[0];
sx q[0];
rz(-2.7543289) q[0];
sx q[0];
rz(1.467147) q[0];
rz(-1.4779429) q[1];
sx q[1];
rz(-1.7324305) q[1];
sx q[1];
rz(-1.0999701) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3008793) q[0];
sx q[0];
rz(-1.3726808) q[0];
sx q[0];
rz(-1.677657) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2312382) q[2];
sx q[2];
rz(-1.8446088) q[2];
sx q[2];
rz(1.6664315) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8588455) q[1];
sx q[1];
rz(-1.2371776) q[1];
sx q[1];
rz(-0.40550123) q[1];
rz(-pi) q[2];
rz(2.1446225) q[3];
sx q[3];
rz(-1.1862635) q[3];
sx q[3];
rz(-1.1763193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.786342) q[2];
sx q[2];
rz(-2.2189271) q[2];
sx q[2];
rz(1.5938909) q[2];
rz(-1.330438) q[3];
sx q[3];
rz(-1.3733613) q[3];
sx q[3];
rz(2.0681341) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.380577) q[0];
sx q[0];
rz(-1.8208193) q[0];
sx q[0];
rz(0.15783489) q[0];
rz(0.24179587) q[1];
sx q[1];
rz(-0.38547412) q[1];
sx q[1];
rz(-2.6604624) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4241971) q[0];
sx q[0];
rz(-1.4823874) q[0];
sx q[0];
rz(0.7128977) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.35901661) q[2];
sx q[2];
rz(-1.8098157) q[2];
sx q[2];
rz(2.0109107) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.85192902) q[1];
sx q[1];
rz(-0.73411513) q[1];
sx q[1];
rz(-1.0490225) q[1];
x q[2];
rz(-1.8577544) q[3];
sx q[3];
rz(-1.7718414) q[3];
sx q[3];
rz(1.505681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.87159291) q[2];
sx q[2];
rz(-0.81120482) q[2];
sx q[2];
rz(2.8672186) q[2];
rz(1.4659878) q[3];
sx q[3];
rz(-1.2803187) q[3];
sx q[3];
rz(-0.93332851) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45193732) q[0];
sx q[0];
rz(-0.87757293) q[0];
sx q[0];
rz(1.0614606) q[0];
rz(2.6929216) q[1];
sx q[1];
rz(-0.45495382) q[1];
sx q[1];
rz(-1.7015069) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7913623) q[0];
sx q[0];
rz(-1.4141948) q[0];
sx q[0];
rz(-0.76275218) q[0];
x q[1];
rz(-1.1665197) q[2];
sx q[2];
rz(-1.4540404) q[2];
sx q[2];
rz(-1.5424042) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6430407) q[1];
sx q[1];
rz(-1.3128377) q[1];
sx q[1];
rz(-0.027679701) q[1];
x q[2];
rz(0.10847059) q[3];
sx q[3];
rz(-2.1785695) q[3];
sx q[3];
rz(-2.8911107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3638641) q[2];
sx q[2];
rz(-1.8455467) q[2];
sx q[2];
rz(1.1931922) q[2];
rz(0.39919546) q[3];
sx q[3];
rz(-1.4123071) q[3];
sx q[3];
rz(0.027755888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4611918) q[0];
sx q[0];
rz(-1.2874648) q[0];
sx q[0];
rz(0.68049085) q[0];
rz(2.3091799) q[1];
sx q[1];
rz(-1.8871555) q[1];
sx q[1];
rz(0.15636538) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.098861) q[0];
sx q[0];
rz(-1.2953238) q[0];
sx q[0];
rz(-1.5163438) q[0];
rz(0.6758718) q[2];
sx q[2];
rz(-2.1891433) q[2];
sx q[2];
rz(-0.72406521) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.63850832) q[1];
sx q[1];
rz(-1.2933532) q[1];
sx q[1];
rz(1.2289407) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0632319) q[3];
sx q[3];
rz(-0.2476736) q[3];
sx q[3];
rz(1.7223235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.51123315) q[2];
sx q[2];
rz(-2.4801621) q[2];
sx q[2];
rz(-2.2516001) q[2];
rz(-2.6651799) q[3];
sx q[3];
rz(-2.2178631) q[3];
sx q[3];
rz(-2.711003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42590672) q[0];
sx q[0];
rz(-1.3193193) q[0];
sx q[0];
rz(-2.529378) q[0];
rz(1.3587492) q[1];
sx q[1];
rz(-2.3814059) q[1];
sx q[1];
rz(-2.8403958) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4801203) q[0];
sx q[0];
rz(-1.299322) q[0];
sx q[0];
rz(-0.090711509) q[0];
rz(-pi) q[1];
rz(1.1015755) q[2];
sx q[2];
rz(-2.42832) q[2];
sx q[2];
rz(0.81339806) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0750904) q[1];
sx q[1];
rz(-1.9478223) q[1];
sx q[1];
rz(2.6860363) q[1];
x q[2];
rz(0.73681085) q[3];
sx q[3];
rz(-2.7932146) q[3];
sx q[3];
rz(-2.762037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.23413868) q[2];
sx q[2];
rz(-2.335304) q[2];
sx q[2];
rz(-1.0478728) q[2];
rz(-2.786484) q[3];
sx q[3];
rz(-0.57098782) q[3];
sx q[3];
rz(2.4560438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.909914) q[0];
sx q[0];
rz(-0.34220085) q[0];
sx q[0];
rz(-2.8177596) q[0];
rz(0.58770761) q[1];
sx q[1];
rz(-0.67981845) q[1];
sx q[1];
rz(0.30276611) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97070951) q[0];
sx q[0];
rz(-1.8237178) q[0];
sx q[0];
rz(2.80632) q[0];
x q[1];
rz(0.23719533) q[2];
sx q[2];
rz(-1.5233808) q[2];
sx q[2];
rz(-0.42979017) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0787133) q[1];
sx q[1];
rz(-0.67725607) q[1];
sx q[1];
rz(0.96911624) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0400502) q[3];
sx q[3];
rz(-2.663739) q[3];
sx q[3];
rz(0.32020204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4964464) q[2];
sx q[2];
rz(-2.1108997) q[2];
sx q[2];
rz(2.5808064) q[2];
rz(1.0049817) q[3];
sx q[3];
rz(-1.2882261) q[3];
sx q[3];
rz(1.0952449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.0305369) q[0];
sx q[0];
rz(-2.472214) q[0];
sx q[0];
rz(-2.3916767) q[0];
rz(-0.38052446) q[1];
sx q[1];
rz(-2.1546202) q[1];
sx q[1];
rz(-2.1662625) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.014054) q[0];
sx q[0];
rz(-0.97146266) q[0];
sx q[0];
rz(2.4571193) q[0];
rz(1.3644977) q[2];
sx q[2];
rz(-2.8562284) q[2];
sx q[2];
rz(2.9138164) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.9561466) q[1];
sx q[1];
rz(-0.6022033) q[1];
sx q[1];
rz(0.57149296) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1998243) q[3];
sx q[3];
rz(-0.36702752) q[3];
sx q[3];
rz(-1.1197156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9425977) q[2];
sx q[2];
rz(-2.2364605) q[2];
sx q[2];
rz(1.0814166) q[2];
rz(0.069247581) q[3];
sx q[3];
rz(-1.705575) q[3];
sx q[3];
rz(-0.92250219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0104495) q[0];
sx q[0];
rz(-2.8157225) q[0];
sx q[0];
rz(-0.3057873) q[0];
rz(0.30793134) q[1];
sx q[1];
rz(-1.6653857) q[1];
sx q[1];
rz(1.9889132) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4534396) q[0];
sx q[0];
rz(-1.3092293) q[0];
sx q[0];
rz(0.75193172) q[0];
x q[1];
rz(0.20736097) q[2];
sx q[2];
rz(-0.93610686) q[2];
sx q[2];
rz(0.35516741) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6438442) q[1];
sx q[1];
rz(-1.1670928) q[1];
sx q[1];
rz(0.35404842) q[1];
rz(-1.5074499) q[3];
sx q[3];
rz(-1.5858486) q[3];
sx q[3];
rz(1.2069595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.61047381) q[2];
sx q[2];
rz(-1.0627154) q[2];
sx q[2];
rz(-1.6884241) q[2];
rz(-2.6609663) q[3];
sx q[3];
rz(-2.5745945) q[3];
sx q[3];
rz(1.7467197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65344812) q[0];
sx q[0];
rz(-1.6076037) q[0];
sx q[0];
rz(-1.6758767) q[0];
rz(-2.0987971) q[1];
sx q[1];
rz(-1.4029618) q[1];
sx q[1];
rz(-0.89697368) q[1];
rz(-2.7077012) q[2];
sx q[2];
rz(-0.57882751) q[2];
sx q[2];
rz(-0.71633518) q[2];
rz(1.5299464) q[3];
sx q[3];
rz(-1.2129285) q[3];
sx q[3];
rz(-0.34679888) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
