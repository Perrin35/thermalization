OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.3192531) q[0];
sx q[0];
rz(-2.1847794) q[0];
sx q[0];
rz(-1.4332888) q[0];
rz(-0.23694555) q[1];
sx q[1];
rz(-1.911093) q[1];
sx q[1];
rz(1.9967611) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0639609) q[0];
sx q[0];
rz(-1.6368757) q[0];
sx q[0];
rz(1.328701) q[0];
rz(-pi) q[1];
rz(1.1510017) q[2];
sx q[2];
rz(-1.3090054) q[2];
sx q[2];
rz(-2.7671438) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.45935985) q[1];
sx q[1];
rz(-2.0233005) q[1];
sx q[1];
rz(-3.1205936) q[1];
x q[2];
rz(1.8688698) q[3];
sx q[3];
rz(-1.0469336) q[3];
sx q[3];
rz(-3.1071172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0248489) q[2];
sx q[2];
rz(-1.3329093) q[2];
sx q[2];
rz(-1.2343181) q[2];
rz(-1.0553137) q[3];
sx q[3];
rz(-1.0588667) q[3];
sx q[3];
rz(1.9887259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9556483) q[0];
sx q[0];
rz(-1.9444822) q[0];
sx q[0];
rz(2.3117075) q[0];
rz(-0.25289598) q[1];
sx q[1];
rz(-0.77650944) q[1];
sx q[1];
rz(-2.2944962) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23628274) q[0];
sx q[0];
rz(-0.61203996) q[0];
sx q[0];
rz(0.73487868) q[0];
x q[1];
rz(-0.29909889) q[2];
sx q[2];
rz(-0.80183376) q[2];
sx q[2];
rz(-2.9849844) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8250679) q[1];
sx q[1];
rz(-1.0267338) q[1];
sx q[1];
rz(-1.6506509) q[1];
x q[2];
rz(-0.91244016) q[3];
sx q[3];
rz(-2.6077301) q[3];
sx q[3];
rz(-2.562059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.117924) q[2];
sx q[2];
rz(-2.7682722) q[2];
sx q[2];
rz(-2.4222597) q[2];
rz(1.8524648) q[3];
sx q[3];
rz(-0.23580655) q[3];
sx q[3];
rz(1.4891362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8711202) q[0];
sx q[0];
rz(-1.2501165) q[0];
sx q[0];
rz(0.57139325) q[0];
rz(-2.1748623) q[1];
sx q[1];
rz(-1.2582018) q[1];
sx q[1];
rz(2.6142696) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2505328) q[0];
sx q[0];
rz(-2.8322304) q[0];
sx q[0];
rz(-3.1191349) q[0];
x q[1];
rz(1.2450561) q[2];
sx q[2];
rz(-0.92391652) q[2];
sx q[2];
rz(-1.5145472) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5273683) q[1];
sx q[1];
rz(-0.80263153) q[1];
sx q[1];
rz(-2.7781092) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8970044) q[3];
sx q[3];
rz(-0.86498125) q[3];
sx q[3];
rz(-1.0969485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.30806914) q[2];
sx q[2];
rz(-1.6922733) q[2];
sx q[2];
rz(0.72675881) q[2];
rz(0.99749342) q[3];
sx q[3];
rz(-0.62524978) q[3];
sx q[3];
rz(1.5184901) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8961287) q[0];
sx q[0];
rz(-1.6587057) q[0];
sx q[0];
rz(-2.1380651) q[0];
rz(-3.1009122) q[1];
sx q[1];
rz(-2.0122806) q[1];
sx q[1];
rz(2.3153268) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6711133) q[0];
sx q[0];
rz(-0.90958909) q[0];
sx q[0];
rz(-0.83755042) q[0];
x q[1];
rz(0.70010186) q[2];
sx q[2];
rz(-1.6119909) q[2];
sx q[2];
rz(0.61566478) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.5845881) q[1];
sx q[1];
rz(-1.6465721) q[1];
sx q[1];
rz(-1.2332548) q[1];
rz(1.1816979) q[3];
sx q[3];
rz(-2.1660921) q[3];
sx q[3];
rz(-0.94483313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5740009) q[2];
sx q[2];
rz(-1.2446128) q[2];
sx q[2];
rz(2.2272002) q[2];
rz(3.0363723) q[3];
sx q[3];
rz(-1.6243694) q[3];
sx q[3];
rz(-0.58925327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88702622) q[0];
sx q[0];
rz(-1.6396602) q[0];
sx q[0];
rz(1.1337093) q[0];
rz(0.0056313593) q[1];
sx q[1];
rz(-1.7040323) q[1];
sx q[1];
rz(-2.5240135) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.782398) q[0];
sx q[0];
rz(-3.087128) q[0];
sx q[0];
rz(0.99476238) q[0];
rz(-pi) q[1];
rz(1.1281625) q[2];
sx q[2];
rz(-2.1652522) q[2];
sx q[2];
rz(-0.15304676) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4434112) q[1];
sx q[1];
rz(-0.65084208) q[1];
sx q[1];
rz(2.0425668) q[1];
rz(3.0693552) q[3];
sx q[3];
rz(-2.9655955) q[3];
sx q[3];
rz(-0.76398677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.79409838) q[2];
sx q[2];
rz(-1.2604875) q[2];
sx q[2];
rz(2.9079672) q[2];
rz(-0.99308333) q[3];
sx q[3];
rz(-0.94998327) q[3];
sx q[3];
rz(-0.63794199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96930209) q[0];
sx q[0];
rz(-3.0811221) q[0];
sx q[0];
rz(-0.0897952) q[0];
rz(0.88109294) q[1];
sx q[1];
rz(-2.6125364) q[1];
sx q[1];
rz(2.9615013) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20443944) q[0];
sx q[0];
rz(-0.75980543) q[0];
sx q[0];
rz(-0.8888437) q[0];
x q[1];
rz(-2.9526268) q[2];
sx q[2];
rz(-0.97550387) q[2];
sx q[2];
rz(1.8408066) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.99299586) q[1];
sx q[1];
rz(-2.8088048) q[1];
sx q[1];
rz(-0.51698835) q[1];
rz(1.5363541) q[3];
sx q[3];
rz(-0.4930217) q[3];
sx q[3];
rz(2.8480414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.883541) q[2];
sx q[2];
rz(-1.5641314) q[2];
sx q[2];
rz(1.9656666) q[2];
rz(-3.0684493) q[3];
sx q[3];
rz(-0.72435838) q[3];
sx q[3];
rz(-0.49889645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95773762) q[0];
sx q[0];
rz(-2.4860005) q[0];
sx q[0];
rz(1.7234329) q[0];
rz(1.9765967) q[1];
sx q[1];
rz(-0.55924758) q[1];
sx q[1];
rz(-0.040239008) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2049853) q[0];
sx q[0];
rz(-1.6523223) q[0];
sx q[0];
rz(1.3315014) q[0];
x q[1];
rz(-1.1646541) q[2];
sx q[2];
rz(-1.6132266) q[2];
sx q[2];
rz(0.050780642) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.38899598) q[1];
sx q[1];
rz(-0.95655671) q[1];
sx q[1];
rz(-0.1715626) q[1];
x q[2];
rz(-1.8764898) q[3];
sx q[3];
rz(-1.2672658) q[3];
sx q[3];
rz(-0.85206735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.13850257) q[2];
sx q[2];
rz(-2.3562727) q[2];
sx q[2];
rz(2.9677532) q[2];
rz(1.396817) q[3];
sx q[3];
rz(-1.4649748) q[3];
sx q[3];
rz(-0.85913908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.023507) q[0];
sx q[0];
rz(-0.21723391) q[0];
sx q[0];
rz(-1.404495) q[0];
rz(-1.9346168) q[1];
sx q[1];
rz(-1.4852306) q[1];
sx q[1];
rz(1.5054024) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0818427) q[0];
sx q[0];
rz(-1.8847701) q[0];
sx q[0];
rz(0.44072515) q[0];
rz(-2.2340441) q[2];
sx q[2];
rz(-2.2758256) q[2];
sx q[2];
rz(-2.7570587) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1464403) q[1];
sx q[1];
rz(-2.4447828) q[1];
sx q[1];
rz(-0.050683024) q[1];
rz(-pi) q[2];
rz(-2.8137915) q[3];
sx q[3];
rz(-0.9947239) q[3];
sx q[3];
rz(-0.93186659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.76190844) q[2];
sx q[2];
rz(-0.48019335) q[2];
sx q[2];
rz(2.9910679) q[2];
rz(-1.5395509) q[3];
sx q[3];
rz(-1.1952885) q[3];
sx q[3];
rz(-2.5185744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4432916) q[0];
sx q[0];
rz(-1.1346096) q[0];
sx q[0];
rz(2.495893) q[0];
rz(2.6422016) q[1];
sx q[1];
rz(-1.4285409) q[1];
sx q[1];
rz(1.8766778) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18683534) q[0];
sx q[0];
rz(-1.1823913) q[0];
sx q[0];
rz(-0.79083058) q[0];
x q[1];
rz(2.0090946) q[2];
sx q[2];
rz(-2.0813137) q[2];
sx q[2];
rz(2.8709708) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7281108) q[1];
sx q[1];
rz(-2.9915447) q[1];
sx q[1];
rz(0.24197443) q[1];
rz(-pi) q[2];
rz(2.7183652) q[3];
sx q[3];
rz(-1.3864577) q[3];
sx q[3];
rz(-1.1545899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.50679961) q[2];
sx q[2];
rz(-2.1719077) q[2];
sx q[2];
rz(0.52465049) q[2];
rz(-2.5850885) q[3];
sx q[3];
rz(-1.5126712) q[3];
sx q[3];
rz(-2.7867253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5234579) q[0];
sx q[0];
rz(-0.83704656) q[0];
sx q[0];
rz(2.8961704) q[0];
rz(0.55631176) q[1];
sx q[1];
rz(-1.6106771) q[1];
sx q[1];
rz(1.4642749) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7510371) q[0];
sx q[0];
rz(-1.7883736) q[0];
sx q[0];
rz(1.4575973) q[0];
rz(0.22428959) q[2];
sx q[2];
rz(-2.1456492) q[2];
sx q[2];
rz(1.526051) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7732685) q[1];
sx q[1];
rz(-1.8414294) q[1];
sx q[1];
rz(0.86398217) q[1];
rz(0.085300307) q[3];
sx q[3];
rz(-1.7461667) q[3];
sx q[3];
rz(0.264197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2422553) q[2];
sx q[2];
rz(-1.9031886) q[2];
sx q[2];
rz(-0.74550068) q[2];
rz(1.7177104) q[3];
sx q[3];
rz(-2.8427377) q[3];
sx q[3];
rz(-1.1283114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.8650919) q[0];
sx q[0];
rz(-1.6456974) q[0];
sx q[0];
rz(0.67281848) q[0];
rz(-1.2561692) q[1];
sx q[1];
rz(-2.3354463) q[1];
sx q[1];
rz(-1.068402) q[1];
rz(-2.4297498) q[2];
sx q[2];
rz(-1.9528452) q[2];
sx q[2];
rz(-0.14609329) q[2];
rz(-0.57337702) q[3];
sx q[3];
rz(-1.2497414) q[3];
sx q[3];
rz(-2.9336815) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
