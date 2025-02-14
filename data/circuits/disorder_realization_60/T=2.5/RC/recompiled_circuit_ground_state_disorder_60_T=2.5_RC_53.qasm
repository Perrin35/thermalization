OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.085091703) q[0];
sx q[0];
rz(-0.29274517) q[0];
sx q[0];
rz(2.2226287) q[0];
rz(2.7594944) q[1];
sx q[1];
rz(-0.16799071) q[1];
sx q[1];
rz(-1.1408495) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5802683) q[0];
sx q[0];
rz(-1.7467357) q[0];
sx q[0];
rz(-0.12138155) q[0];
x q[1];
rz(0.72854002) q[2];
sx q[2];
rz(-2.72914) q[2];
sx q[2];
rz(-2.0714687) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.029143229) q[1];
sx q[1];
rz(-1.1894798) q[1];
sx q[1];
rz(1.4277532) q[1];
x q[2];
rz(1.857418) q[3];
sx q[3];
rz(-0.45551963) q[3];
sx q[3];
rz(-2.974433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4484278) q[2];
sx q[2];
rz(-2.0824671) q[2];
sx q[2];
rz(-1.1761752) q[2];
rz(-0.042393427) q[3];
sx q[3];
rz(-1.8040801) q[3];
sx q[3];
rz(-2.6068408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46928826) q[0];
sx q[0];
rz(-0.35457087) q[0];
sx q[0];
rz(-0.18445036) q[0];
rz(0.09672673) q[1];
sx q[1];
rz(-1.5238785) q[1];
sx q[1];
rz(1.2299889) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0320659) q[0];
sx q[0];
rz(-1.7735039) q[0];
sx q[0];
rz(2.7687182) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.619009) q[2];
sx q[2];
rz(-2.1297751) q[2];
sx q[2];
rz(0.79254442) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0945541) q[1];
sx q[1];
rz(-0.31289214) q[1];
sx q[1];
rz(-2.439586) q[1];
x q[2];
rz(1.473068) q[3];
sx q[3];
rz(-0.37993452) q[3];
sx q[3];
rz(-1.987628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.35473287) q[2];
sx q[2];
rz(-2.730098) q[2];
sx q[2];
rz(-0.011431781) q[2];
rz(1.8276021) q[3];
sx q[3];
rz(-1.7766137) q[3];
sx q[3];
rz(2.438681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74362022) q[0];
sx q[0];
rz(-3.008606) q[0];
sx q[0];
rz(2.5308894) q[0];
rz(-0.01344219) q[1];
sx q[1];
rz(-1.2862658) q[1];
sx q[1];
rz(0.59649831) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4152348) q[0];
sx q[0];
rz(-1.7371375) q[0];
sx q[0];
rz(-2.1852058) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1488068) q[2];
sx q[2];
rz(-1.3137378) q[2];
sx q[2];
rz(1.3959194) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8104659) q[1];
sx q[1];
rz(-2.2743723) q[1];
sx q[1];
rz(1.1570279) q[1];
x q[2];
rz(3.1085988) q[3];
sx q[3];
rz(-2.1244496) q[3];
sx q[3];
rz(0.43137303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.67768031) q[2];
sx q[2];
rz(-1.8296506) q[2];
sx q[2];
rz(-3.1413191) q[2];
rz(-1.6655946) q[3];
sx q[3];
rz(-2.4725584) q[3];
sx q[3];
rz(-3.0345548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6855455) q[0];
sx q[0];
rz(-2.9750329) q[0];
sx q[0];
rz(1.9786932) q[0];
rz(2.1641425) q[1];
sx q[1];
rz(-1.5403427) q[1];
sx q[1];
rz(-1.5862484) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.259183) q[0];
sx q[0];
rz(-1.5611739) q[0];
sx q[0];
rz(0.0077879328) q[0];
rz(-pi) q[1];
rz(0.78537099) q[2];
sx q[2];
rz(-1.6024688) q[2];
sx q[2];
rz(-2.9065064) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4443497) q[1];
sx q[1];
rz(-1.0740136) q[1];
sx q[1];
rz(1.5543943) q[1];
x q[2];
rz(1.1337676) q[3];
sx q[3];
rz(-2.6053782) q[3];
sx q[3];
rz(1.0325583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.9419452) q[2];
sx q[2];
rz(-0.5117828) q[2];
sx q[2];
rz(-2.9072348) q[2];
rz(2.6394081) q[3];
sx q[3];
rz(-2.1000704) q[3];
sx q[3];
rz(2.485399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.31768826) q[0];
sx q[0];
rz(-2.3264139) q[0];
sx q[0];
rz(1.3932047) q[0];
rz(2.4488917) q[1];
sx q[1];
rz(-1.0857948) q[1];
sx q[1];
rz(-1.0903953) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2176184) q[0];
sx q[0];
rz(-2.0912785) q[0];
sx q[0];
rz(-2.897516) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.52509016) q[2];
sx q[2];
rz(-1.0193362) q[2];
sx q[2];
rz(2.9208825) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.81933669) q[1];
sx q[1];
rz(-2.6988479) q[1];
sx q[1];
rz(-1.2135452) q[1];
rz(-1.8282918) q[3];
sx q[3];
rz(-0.55488834) q[3];
sx q[3];
rz(-2.1438847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.580487) q[2];
sx q[2];
rz(-1.9876391) q[2];
sx q[2];
rz(0.53606501) q[2];
rz(2.8481893) q[3];
sx q[3];
rz(-1.9304201) q[3];
sx q[3];
rz(2.6611879) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90477657) q[0];
sx q[0];
rz(-2.1892956) q[0];
sx q[0];
rz(-3.0425518) q[0];
rz(1.0320484) q[1];
sx q[1];
rz(-2.2910304) q[1];
sx q[1];
rz(1.438407) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0625298) q[0];
sx q[0];
rz(-0.89508104) q[0];
sx q[0];
rz(-0.32846256) q[0];
rz(-0.24865227) q[2];
sx q[2];
rz(-0.2919582) q[2];
sx q[2];
rz(-2.0338634) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.21675303) q[1];
sx q[1];
rz(-2.3509563) q[1];
sx q[1];
rz(-3.1086712) q[1];
x q[2];
rz(-1.4505054) q[3];
sx q[3];
rz(-1.866739) q[3];
sx q[3];
rz(0.048792865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9859887) q[2];
sx q[2];
rz(-0.5039379) q[2];
sx q[2];
rz(0.8209374) q[2];
rz(1.4486897) q[3];
sx q[3];
rz(-2.4235348) q[3];
sx q[3];
rz(1.4887571) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2445225) q[0];
sx q[0];
rz(-0.91461602) q[0];
sx q[0];
rz(0.50265092) q[0];
rz(-1.9082759) q[1];
sx q[1];
rz(-0.93524593) q[1];
sx q[1];
rz(-0.9800235) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4122881) q[0];
sx q[0];
rz(-1.8192023) q[0];
sx q[0];
rz(-0.24459837) q[0];
rz(-pi) q[1];
rz(0.77620164) q[2];
sx q[2];
rz(-1.7817678) q[2];
sx q[2];
rz(-0.38061505) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.087638559) q[1];
sx q[1];
rz(-2.3628937) q[1];
sx q[1];
rz(-2.8092572) q[1];
x q[2];
rz(1.3408889) q[3];
sx q[3];
rz(-2.3977444) q[3];
sx q[3];
rz(-2.6514253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9194453) q[2];
sx q[2];
rz(-1.8541226) q[2];
sx q[2];
rz(1.7670828) q[2];
rz(1.3162656) q[3];
sx q[3];
rz(-0.85597435) q[3];
sx q[3];
rz(2.159806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.027243622) q[0];
sx q[0];
rz(-0.39334941) q[0];
sx q[0];
rz(-2.223176) q[0];
rz(-0.20485993) q[1];
sx q[1];
rz(-1.058895) q[1];
sx q[1];
rz(-1.6580261) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.525104) q[0];
sx q[0];
rz(-1.8601928) q[0];
sx q[0];
rz(0.15286907) q[0];
rz(-pi) q[1];
rz(-0.81088938) q[2];
sx q[2];
rz(-1.4810908) q[2];
sx q[2];
rz(-0.16505884) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.073347884) q[1];
sx q[1];
rz(-2.4363344) q[1];
sx q[1];
rz(1.4346992) q[1];
rz(0.71916716) q[3];
sx q[3];
rz(-1.5336108) q[3];
sx q[3];
rz(2.3553762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.001699) q[2];
sx q[2];
rz(-0.49109444) q[2];
sx q[2];
rz(1.8074544) q[2];
rz(-0.88207465) q[3];
sx q[3];
rz(-2.4460654) q[3];
sx q[3];
rz(-1.2923366) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0598711) q[0];
sx q[0];
rz(-0.43161714) q[0];
sx q[0];
rz(-0.01817848) q[0];
rz(0.40953088) q[1];
sx q[1];
rz(-1.4298226) q[1];
sx q[1];
rz(0.10083625) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6438599) q[0];
sx q[0];
rz(-2.5063305) q[0];
sx q[0];
rz(-3.0833901) q[0];
x q[1];
rz(-0.42068215) q[2];
sx q[2];
rz(-1.6170338) q[2];
sx q[2];
rz(-1.7939523) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.455115) q[1];
sx q[1];
rz(-1.9338738) q[1];
sx q[1];
rz(2.8836714) q[1];
rz(-pi) q[2];
rz(3.0206817) q[3];
sx q[3];
rz(-1.3252581) q[3];
sx q[3];
rz(2.3335378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.33710256) q[2];
sx q[2];
rz(-0.10919315) q[2];
sx q[2];
rz(1.8936554) q[2];
rz(-1.9167871) q[3];
sx q[3];
rz(-2.2962511) q[3];
sx q[3];
rz(-0.065751806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2357904) q[0];
sx q[0];
rz(-1.4757272) q[0];
sx q[0];
rz(-2.8210848) q[0];
rz(0.39101741) q[1];
sx q[1];
rz(-1.1415488) q[1];
sx q[1];
rz(1.549622) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3282916) q[0];
sx q[0];
rz(-1.7696972) q[0];
sx q[0];
rz(1.7651674) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6274849) q[2];
sx q[2];
rz(-1.7127345) q[2];
sx q[2];
rz(-2.1588003) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.95028535) q[1];
sx q[1];
rz(-2.7627594) q[1];
sx q[1];
rz(1.8089507) q[1];
rz(1.5632838) q[3];
sx q[3];
rz(-1.148461) q[3];
sx q[3];
rz(1.0653121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0512507) q[2];
sx q[2];
rz(-2.0917459) q[2];
sx q[2];
rz(0.40038294) q[2];
rz(-2.9054437) q[3];
sx q[3];
rz(-0.103424) q[3];
sx q[3];
rz(-1.0108488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.878933) q[0];
sx q[0];
rz(-2.6317609) q[0];
sx q[0];
rz(1.2608933) q[0];
rz(-0.48514584) q[1];
sx q[1];
rz(-0.79882516) q[1];
sx q[1];
rz(-0.47732236) q[1];
rz(-0.12564364) q[2];
sx q[2];
rz(-1.7775272) q[2];
sx q[2];
rz(0.20128332) q[2];
rz(-1.7988206) q[3];
sx q[3];
rz(-1.6327471) q[3];
sx q[3];
rz(-3.0625797) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
