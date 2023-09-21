OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0302304) q[0];
sx q[0];
rz(5.6279866) q[0];
sx q[0];
rz(10.401166) q[0];
rz(-0.916565) q[1];
sx q[1];
rz(2.38148) q[1];
sx q[1];
rz(10.217477) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7615258) q[0];
sx q[0];
rz(-0.67136429) q[0];
sx q[0];
rz(-1.6094366) q[0];
rz(-pi) q[1];
rz(-1.3408469) q[2];
sx q[2];
rz(-1.0385498) q[2];
sx q[2];
rz(-1.6543653) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.76268643) q[1];
sx q[1];
rz(-1.1211809) q[1];
sx q[1];
rz(-0.080888918) q[1];
rz(0.80134942) q[3];
sx q[3];
rz(-2.1270463) q[3];
sx q[3];
rz(-2.706561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4684326) q[2];
sx q[2];
rz(-0.94753733) q[2];
sx q[2];
rz(2.6426278) q[2];
rz(0.4237825) q[3];
sx q[3];
rz(-1.2200004) q[3];
sx q[3];
rz(-0.012714816) q[3];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.760261) q[0];
sx q[0];
rz(-1.6853764) q[0];
sx q[0];
rz(-0.99200845) q[0];
rz(0.17653067) q[1];
sx q[1];
rz(-2.9361528) q[1];
sx q[1];
rz(-0.36546779) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8651543) q[0];
sx q[0];
rz(-1.3197474) q[0];
sx q[0];
rz(2.6692078) q[0];
x q[1];
rz(2.2720049) q[2];
sx q[2];
rz(-2.3361932) q[2];
sx q[2];
rz(-0.7676917) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5181527) q[1];
sx q[1];
rz(-1.232051) q[1];
sx q[1];
rz(-2.5046196) q[1];
rz(2.9404853) q[3];
sx q[3];
rz(-1.4521986) q[3];
sx q[3];
rz(3.0191641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3084597) q[2];
sx q[2];
rz(-2.5434727) q[2];
sx q[2];
rz(2.0642521) q[2];
rz(0.16263738) q[3];
sx q[3];
rz(-1.2626303) q[3];
sx q[3];
rz(1.8796273) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66115528) q[0];
sx q[0];
rz(-0.85395956) q[0];
sx q[0];
rz(0.61082947) q[0];
rz(-1.4340596) q[1];
sx q[1];
rz(-1.779498) q[1];
sx q[1];
rz(-2.6909713) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.492309) q[0];
sx q[0];
rz(-1.603924) q[0];
sx q[0];
rz(-1.6487728) q[0];
x q[1];
rz(1.9446104) q[2];
sx q[2];
rz(-1.2842368) q[2];
sx q[2];
rz(-2.8754004) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1994985) q[1];
sx q[1];
rz(-2.4220785) q[1];
sx q[1];
rz(2.662563) q[1];
rz(-pi) q[2];
rz(-2.2654128) q[3];
sx q[3];
rz(-1.8091222) q[3];
sx q[3];
rz(0.5660457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0314363) q[2];
sx q[2];
rz(-2.0560196) q[2];
sx q[2];
rz(-1.191656) q[2];
rz(-0.84447652) q[3];
sx q[3];
rz(-0.28544393) q[3];
sx q[3];
rz(-2.5604131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0464756) q[0];
sx q[0];
rz(-1.7549544) q[0];
sx q[0];
rz(-1.1412096) q[0];
rz(-1.8449239) q[1];
sx q[1];
rz(-2.7099687) q[1];
sx q[1];
rz(-2.8362714) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6396219) q[0];
sx q[0];
rz(-1.0143711) q[0];
sx q[0];
rz(-1.200838) q[0];
rz(0.73037578) q[2];
sx q[2];
rz(-0.59130284) q[2];
sx q[2];
rz(3.0808466) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7118453) q[1];
sx q[1];
rz(-1.4652805) q[1];
sx q[1];
rz(2.3464416) q[1];
x q[2];
rz(-0.51586113) q[3];
sx q[3];
rz(-2.0277884) q[3];
sx q[3];
rz(0.36992321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9332463) q[2];
sx q[2];
rz(-1.5657921) q[2];
sx q[2];
rz(-0.34710458) q[2];
rz(-2.7582205) q[3];
sx q[3];
rz(-2.255286) q[3];
sx q[3];
rz(-2.0736407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7326996) q[0];
sx q[0];
rz(-2.2786031) q[0];
sx q[0];
rz(-2.6026759) q[0];
rz(1.7319038) q[1];
sx q[1];
rz(-2.6795645) q[1];
sx q[1];
rz(2.7358823) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.021279871) q[0];
sx q[0];
rz(-0.73182523) q[0];
sx q[0];
rz(-2.7243607) q[0];
rz(-1.1666127) q[2];
sx q[2];
rz(-2.5378072) q[2];
sx q[2];
rz(-0.052978901) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0438784) q[1];
sx q[1];
rz(-1.8227302) q[1];
sx q[1];
rz(1.21752) q[1];
rz(-pi) q[2];
x q[2];
rz(0.76463502) q[3];
sx q[3];
rz(-1.4009762) q[3];
sx q[3];
rz(-3.0499383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6695909) q[2];
sx q[2];
rz(-0.80299032) q[2];
sx q[2];
rz(-1.6873138) q[2];
rz(2.3029095) q[3];
sx q[3];
rz(-1.8554153) q[3];
sx q[3];
rz(-1.7680761) q[3];
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
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31081653) q[0];
sx q[0];
rz(-2.5287479) q[0];
sx q[0];
rz(2.7344761) q[0];
rz(2.282417) q[1];
sx q[1];
rz(-1.5770864) q[1];
sx q[1];
rz(-1.8242594) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.84264) q[0];
sx q[0];
rz(-1.5931804) q[0];
sx q[0];
rz(1.0277261) q[0];
rz(-pi) q[1];
x q[1];
rz(0.98951927) q[2];
sx q[2];
rz(-0.98810722) q[2];
sx q[2];
rz(1.4942102) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.085560872) q[1];
sx q[1];
rz(-0.83562682) q[1];
sx q[1];
rz(-3.1378156) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1727354) q[3];
sx q[3];
rz(-0.67776206) q[3];
sx q[3];
rz(-2.4399151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0417827) q[2];
sx q[2];
rz(-0.87974292) q[2];
sx q[2];
rz(2.2463634) q[2];
rz(1.2081395) q[3];
sx q[3];
rz(-2.4619305) q[3];
sx q[3];
rz(-2.4334548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5773425) q[0];
sx q[0];
rz(-2.0944493) q[0];
sx q[0];
rz(2.0136925) q[0];
rz(-2.9290507) q[1];
sx q[1];
rz(-1.3393341) q[1];
sx q[1];
rz(-1.2316661) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54264442) q[0];
sx q[0];
rz(-2.0751187) q[0];
sx q[0];
rz(1.1619316) q[0];
rz(-pi) q[1];
rz(2.2052843) q[2];
sx q[2];
rz(-2.1116743) q[2];
sx q[2];
rz(2.4782431) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1612138) q[1];
sx q[1];
rz(-1.8002138) q[1];
sx q[1];
rz(-1.4646261) q[1];
rz(1.6024186) q[3];
sx q[3];
rz(-2.5368689) q[3];
sx q[3];
rz(-1.5550169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.82905519) q[2];
sx q[2];
rz(-0.91784358) q[2];
sx q[2];
rz(-0.60116872) q[2];
rz(-0.51856315) q[3];
sx q[3];
rz(-2.8278567) q[3];
sx q[3];
rz(-1.7140088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9940051) q[0];
sx q[0];
rz(-3.1412791) q[0];
sx q[0];
rz(0.073154733) q[0];
rz(1.0271429) q[1];
sx q[1];
rz(-1.606753) q[1];
sx q[1];
rz(2.5491319) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52623977) q[0];
sx q[0];
rz(-2.549299) q[0];
sx q[0];
rz(3.0583529) q[0];
rz(-pi) q[1];
rz(-1.0560889) q[2];
sx q[2];
rz(-2.0851496) q[2];
sx q[2];
rz(2.8714542) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.29443024) q[1];
sx q[1];
rz(-1.8995598) q[1];
sx q[1];
rz(2.2734103) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9755649) q[3];
sx q[3];
rz(-1.824607) q[3];
sx q[3];
rz(-2.7491731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5804194) q[2];
sx q[2];
rz(-1.6757123) q[2];
sx q[2];
rz(2.6994761) q[2];
rz(2.2413065) q[3];
sx q[3];
rz(-2.0787652) q[3];
sx q[3];
rz(-3.1282848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1135547) q[0];
sx q[0];
rz(-0.91911626) q[0];
sx q[0];
rz(-1.8369209) q[0];
rz(1.4470944) q[1];
sx q[1];
rz(-1.9686952) q[1];
sx q[1];
rz(-0.56217271) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8174443) q[0];
sx q[0];
rz(-1.9262909) q[0];
sx q[0];
rz(0.72974156) q[0];
rz(1.7299037) q[2];
sx q[2];
rz(-0.53630398) q[2];
sx q[2];
rz(-1.2250587) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0253898) q[1];
sx q[1];
rz(-0.89239489) q[1];
sx q[1];
rz(-0.62979001) q[1];
x q[2];
rz(0.77769827) q[3];
sx q[3];
rz(-0.54293984) q[3];
sx q[3];
rz(-0.87399769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.28329453) q[2];
sx q[2];
rz(-1.7389745) q[2];
sx q[2];
rz(2.1600058) q[2];
rz(0.0043407241) q[3];
sx q[3];
rz(-2.0948295) q[3];
sx q[3];
rz(2.5481352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3146064) q[0];
sx q[0];
rz(-1.2207323) q[0];
sx q[0];
rz(-3.0944968) q[0];
rz(-1.754952) q[1];
sx q[1];
rz(-1.9853233) q[1];
sx q[1];
rz(-0.047853619) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76059607) q[0];
sx q[0];
rz(-2.2254235) q[0];
sx q[0];
rz(-0.036548338) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.677703) q[2];
sx q[2];
rz(-0.44471834) q[2];
sx q[2];
rz(2.5267548) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.25803265) q[1];
sx q[1];
rz(-1.1415392) q[1];
sx q[1];
rz(2.2116823) q[1];
rz(-pi) q[2];
rz(1.5157394) q[3];
sx q[3];
rz(-1.1686472) q[3];
sx q[3];
rz(2.1189342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.64426595) q[2];
sx q[2];
rz(-1.9571807) q[2];
sx q[2];
rz(2.2052374) q[2];
rz(2.3406773) q[3];
sx q[3];
rz(-2.622486) q[3];
sx q[3];
rz(0.67805964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3180278) q[0];
sx q[0];
rz(-0.93946539) q[0];
sx q[0];
rz(-2.9325624) q[0];
rz(0.50280747) q[1];
sx q[1];
rz(-1.8223169) q[1];
sx q[1];
rz(1.0261789) q[1];
rz(0.035632523) q[2];
sx q[2];
rz(-2.590066) q[2];
sx q[2];
rz(-1.8077015) q[2];
rz(-1.1669284) q[3];
sx q[3];
rz(-1.9098017) q[3];
sx q[3];
rz(-3.1168934) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];