OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5366323) q[0];
sx q[0];
rz(-0.16083117) q[0];
sx q[0];
rz(-0.39128006) q[0];
rz(-0.04303509) q[1];
sx q[1];
rz(4.7624762) q[1];
sx q[1];
rz(12.227992) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94253263) q[0];
sx q[0];
rz(-1.8359658) q[0];
sx q[0];
rz(1.2808179) q[0];
x q[1];
rz(-3.0909371) q[2];
sx q[2];
rz(-1.3565552) q[2];
sx q[2];
rz(-0.07344499) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1115626) q[1];
sx q[1];
rz(-1.7610368) q[1];
sx q[1];
rz(-2.2493786) q[1];
rz(-pi) q[2];
rz(1.4621327) q[3];
sx q[3];
rz(-1.1564643) q[3];
sx q[3];
rz(-1.3760096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7732064) q[2];
sx q[2];
rz(-2.1619004) q[2];
sx q[2];
rz(1.0580019) q[2];
rz(0.10231415) q[3];
sx q[3];
rz(-1.5240069) q[3];
sx q[3];
rz(-1.4003096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31829396) q[0];
sx q[0];
rz(-0.75495356) q[0];
sx q[0];
rz(-2.9276983) q[0];
rz(1.8077883) q[1];
sx q[1];
rz(-2.6413481) q[1];
sx q[1];
rz(0.28876367) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28604668) q[0];
sx q[0];
rz(-0.65546765) q[0];
sx q[0];
rz(1.855887) q[0];
rz(-pi) q[1];
rz(2.4756672) q[2];
sx q[2];
rz(-2.1873432) q[2];
sx q[2];
rz(1.2135972) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3067291) q[1];
sx q[1];
rz(-0.22278654) q[1];
sx q[1];
rz(-1.7053717) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.706224) q[3];
sx q[3];
rz(-1.8791128) q[3];
sx q[3];
rz(-1.8086019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5220962) q[2];
sx q[2];
rz(-2.3213826) q[2];
sx q[2];
rz(-1.2843457) q[2];
rz(1.8340825) q[3];
sx q[3];
rz(-1.8239572) q[3];
sx q[3];
rz(2.4431084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(2.4035325) q[0];
sx q[0];
rz(-1.0665251) q[0];
sx q[0];
rz(-2.2962978) q[0];
rz(-2.2360133) q[1];
sx q[1];
rz(-2.5366492) q[1];
sx q[1];
rz(1.1776498) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.986341) q[0];
sx q[0];
rz(-3.0444859) q[0];
sx q[0];
rz(2.9349021) q[0];
rz(-2.3272456) q[2];
sx q[2];
rz(-2.4422514) q[2];
sx q[2];
rz(1.4171079) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.48124734) q[1];
sx q[1];
rz(-1.3637017) q[1];
sx q[1];
rz(0.58086632) q[1];
rz(-0.92567851) q[3];
sx q[3];
rz(-1.4429907) q[3];
sx q[3];
rz(-0.61193454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5251069) q[2];
sx q[2];
rz(-0.944204) q[2];
sx q[2];
rz(0.13060972) q[2];
rz(-1.7836001) q[3];
sx q[3];
rz(-1.0087174) q[3];
sx q[3];
rz(-1.85359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4865049) q[0];
sx q[0];
rz(-2.8772652) q[0];
sx q[0];
rz(-0.41410145) q[0];
rz(0.86530238) q[1];
sx q[1];
rz(-1.2369786) q[1];
sx q[1];
rz(0.6032595) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7287528) q[0];
sx q[0];
rz(-1.9590634) q[0];
sx q[0];
rz(0.30465841) q[0];
rz(-pi) q[1];
rz(-1.499699) q[2];
sx q[2];
rz(-1.9219766) q[2];
sx q[2];
rz(0.49654135) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.016876) q[1];
sx q[1];
rz(-1.6137329) q[1];
sx q[1];
rz(2.7932568) q[1];
x q[2];
rz(1.8835265) q[3];
sx q[3];
rz(-1.4566753) q[3];
sx q[3];
rz(3.0437226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6174751) q[2];
sx q[2];
rz(-1.6092499) q[2];
sx q[2];
rz(2.4821846) q[2];
rz(-3.009033) q[3];
sx q[3];
rz(-2.5949251) q[3];
sx q[3];
rz(-2.4373655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7947614) q[0];
sx q[0];
rz(-2.1857388) q[0];
sx q[0];
rz(0.26926789) q[0];
rz(1.5003834) q[1];
sx q[1];
rz(-1.2790479) q[1];
sx q[1];
rz(-1.0795116) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74306946) q[0];
sx q[0];
rz(-0.67248225) q[0];
sx q[0];
rz(0.29015707) q[0];
rz(-pi) q[1];
rz(0.90957162) q[2];
sx q[2];
rz(-0.9098297) q[2];
sx q[2];
rz(1.329042) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1126571) q[1];
sx q[1];
rz(-1.9017692) q[1];
sx q[1];
rz(-2.079936) q[1];
x q[2];
rz(2.5629839) q[3];
sx q[3];
rz(-0.37453412) q[3];
sx q[3];
rz(2.7079564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.72011224) q[2];
sx q[2];
rz(-0.58028996) q[2];
sx q[2];
rz(-0.87023467) q[2];
rz(-1.8057711) q[3];
sx q[3];
rz(-1.8060874) q[3];
sx q[3];
rz(-0.67572063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0077591) q[0];
sx q[0];
rz(-3.1133911) q[0];
sx q[0];
rz(2.5303685) q[0];
rz(-2.4080343) q[1];
sx q[1];
rz(-1.6707784) q[1];
sx q[1];
rz(2.7582817) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0867423) q[0];
sx q[0];
rz(-2.2630458) q[0];
sx q[0];
rz(1.9350497) q[0];
rz(-pi) q[1];
rz(-2.9401851) q[2];
sx q[2];
rz(-0.83759356) q[2];
sx q[2];
rz(-0.3647441) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9360507) q[1];
sx q[1];
rz(-0.65434736) q[1];
sx q[1];
rz(-1.3597354) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1556968) q[3];
sx q[3];
rz(-2.2533247) q[3];
sx q[3];
rz(2.0180338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.40206259) q[2];
sx q[2];
rz(-0.5548839) q[2];
sx q[2];
rz(-1.0129207) q[2];
rz(-2.5439475) q[3];
sx q[3];
rz(-1.7326071) q[3];
sx q[3];
rz(-0.92596084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0830072) q[0];
sx q[0];
rz(-0.42917955) q[0];
sx q[0];
rz(-0.035932628) q[0];
rz(-1.2220471) q[1];
sx q[1];
rz(-0.67600328) q[1];
sx q[1];
rz(0.24807182) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1386458) q[0];
sx q[0];
rz(-0.32297501) q[0];
sx q[0];
rz(0.92355152) q[0];
rz(2.8808458) q[2];
sx q[2];
rz(-2.6364821) q[2];
sx q[2];
rz(-1.3253761) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8442432) q[1];
sx q[1];
rz(-1.9378758) q[1];
sx q[1];
rz(0.17508164) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5682427) q[3];
sx q[3];
rz(-1.4658584) q[3];
sx q[3];
rz(1.6014639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0566569) q[2];
sx q[2];
rz(-1.0597022) q[2];
sx q[2];
rz(2.6094931) q[2];
rz(-2.687124) q[3];
sx q[3];
rz(-1.5458509) q[3];
sx q[3];
rz(1.2940297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4249975) q[0];
sx q[0];
rz(-1.0497365) q[0];
sx q[0];
rz(-0.24018921) q[0];
rz(2.5364618) q[1];
sx q[1];
rz(-1.3477707) q[1];
sx q[1];
rz(-2.4952707) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8511348) q[0];
sx q[0];
rz(-0.86800985) q[0];
sx q[0];
rz(2.7443307) q[0];
x q[1];
rz(-1.4736299) q[2];
sx q[2];
rz(-2.7787913) q[2];
sx q[2];
rz(1.3902612) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.60257116) q[1];
sx q[1];
rz(-0.50316167) q[1];
sx q[1];
rz(0.43518592) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2262832) q[3];
sx q[3];
rz(-1.3983677) q[3];
sx q[3];
rz(0.61911914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0295082) q[2];
sx q[2];
rz(-1.4936451) q[2];
sx q[2];
rz(-3.0599111) q[2];
rz(2.2937842) q[3];
sx q[3];
rz(-0.40926465) q[3];
sx q[3];
rz(0.22296396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34143701) q[0];
sx q[0];
rz(-2.2969022) q[0];
sx q[0];
rz(-0.89168125) q[0];
rz(1.910123) q[1];
sx q[1];
rz(-2.6530118) q[1];
sx q[1];
rz(-0.60985342) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0783892) q[0];
sx q[0];
rz(-1.1097621) q[0];
sx q[0];
rz(1.0151063) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.78471176) q[2];
sx q[2];
rz(-1.5605188) q[2];
sx q[2];
rz(2.1267872) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.66850502) q[1];
sx q[1];
rz(-2.3021295) q[1];
sx q[1];
rz(-1.1727612) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5242349) q[3];
sx q[3];
rz(-0.34544975) q[3];
sx q[3];
rz(0.54333401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8651198) q[2];
sx q[2];
rz(-0.96729326) q[2];
sx q[2];
rz(-2.8709732) q[2];
rz(0.5451777) q[3];
sx q[3];
rz(-1.3278278) q[3];
sx q[3];
rz(0.24817185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4973064) q[0];
sx q[0];
rz(-1.8338642) q[0];
sx q[0];
rz(-2.9551031) q[0];
rz(-1.9322152) q[1];
sx q[1];
rz(-1.3554363) q[1];
sx q[1];
rz(0.019651042) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3075382) q[0];
sx q[0];
rz(-2.3702633) q[0];
sx q[0];
rz(-2.176009) q[0];
rz(-1.8397464) q[2];
sx q[2];
rz(-1.3185314) q[2];
sx q[2];
rz(2.4076574) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.84148592) q[1];
sx q[1];
rz(-1.7238614) q[1];
sx q[1];
rz(-0.35565214) q[1];
rz(0.76373617) q[3];
sx q[3];
rz(-1.0794355) q[3];
sx q[3];
rz(-1.9869252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.43202117) q[2];
sx q[2];
rz(-0.95546466) q[2];
sx q[2];
rz(1.5891937) q[2];
rz(-3.0324557) q[3];
sx q[3];
rz(-1.8692632) q[3];
sx q[3];
rz(-0.2763589) q[3];
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
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6266262) q[0];
sx q[0];
rz(-1.4780541) q[0];
sx q[0];
rz(1.9274101) q[0];
rz(-1.7264438) q[1];
sx q[1];
rz(-2.1198004) q[1];
sx q[1];
rz(-1.459495) q[1];
rz(-0.44533163) q[2];
sx q[2];
rz(-1.3117322) q[2];
sx q[2];
rz(-1.7016254) q[2];
rz(2.9045803) q[3];
sx q[3];
rz(-0.3921628) q[3];
sx q[3];
rz(-1.8054448) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
