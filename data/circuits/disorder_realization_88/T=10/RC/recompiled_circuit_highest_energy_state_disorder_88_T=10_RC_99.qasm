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
rz(-0.1470546) q[0];
sx q[0];
rz(-1.4015863) q[0];
sx q[0];
rz(0.92010486) q[0];
rz(3.0609581) q[1];
sx q[1];
rz(-0.57890761) q[1];
sx q[1];
rz(-0.97715598) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1950705) q[0];
sx q[0];
rz(-1.3597288) q[0];
sx q[0];
rz(-0.37977438) q[0];
x q[1];
rz(3.009367) q[2];
sx q[2];
rz(-1.5572845) q[2];
sx q[2];
rz(0.15541645) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5861533) q[1];
sx q[1];
rz(-1.7057014) q[1];
sx q[1];
rz(-2.5389266) q[1];
x q[2];
rz(2.7746088) q[3];
sx q[3];
rz(-2.4189848) q[3];
sx q[3];
rz(-2.5740576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5795634) q[2];
sx q[2];
rz(-1.9271489) q[2];
sx q[2];
rz(-2.5207632) q[2];
rz(-2.6260455) q[3];
sx q[3];
rz(-1.4225057) q[3];
sx q[3];
rz(-0.8425042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(1.2466549) q[0];
sx q[0];
rz(-1.5351013) q[0];
sx q[0];
rz(0.57902336) q[0];
rz(-2.31965) q[1];
sx q[1];
rz(-1.6252981) q[1];
sx q[1];
rz(-2.6878405) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0653719) q[0];
sx q[0];
rz(-0.95673086) q[0];
sx q[0];
rz(1.2943511) q[0];
rz(-2.6751509) q[2];
sx q[2];
rz(-1.6714408) q[2];
sx q[2];
rz(-2.538344) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.380341) q[1];
sx q[1];
rz(-0.56219343) q[1];
sx q[1];
rz(0.23558771) q[1];
rz(-0.445153) q[3];
sx q[3];
rz(-1.0694519) q[3];
sx q[3];
rz(-1.2847163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.73432505) q[2];
sx q[2];
rz(-3.0027323) q[2];
sx q[2];
rz(-2.5767051) q[2];
rz(-2.7867553) q[3];
sx q[3];
rz(-0.9762888) q[3];
sx q[3];
rz(-1.7179276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(1.1693717) q[0];
sx q[0];
rz(-2.9300949) q[0];
sx q[0];
rz(0.55150223) q[0];
rz(-0.36477271) q[1];
sx q[1];
rz(-0.33939895) q[1];
sx q[1];
rz(0.50484467) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6423787) q[0];
sx q[0];
rz(-1.9758245) q[0];
sx q[0];
rz(1.3971055) q[0];
rz(-pi) q[1];
rz(-0.41270035) q[2];
sx q[2];
rz(-1.6566212) q[2];
sx q[2];
rz(0.43666652) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.80124679) q[1];
sx q[1];
rz(-2.1643504) q[1];
sx q[1];
rz(-2.1562063) q[1];
rz(0.48074333) q[3];
sx q[3];
rz(-1.9078443) q[3];
sx q[3];
rz(2.2729006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.43368936) q[2];
sx q[2];
rz(-1.7283551) q[2];
sx q[2];
rz(-0.047253963) q[2];
rz(-0.83682483) q[3];
sx q[3];
rz(-0.72332007) q[3];
sx q[3];
rz(2.1164472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4243917) q[0];
sx q[0];
rz(-1.0294788) q[0];
sx q[0];
rz(1.3264054) q[0];
rz(1.4855509) q[1];
sx q[1];
rz(-1.0842208) q[1];
sx q[1];
rz(2.5066689) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0795006) q[0];
sx q[0];
rz(-1.5720815) q[0];
sx q[0];
rz(1.5888402) q[0];
rz(-pi) q[1];
x q[1];
rz(0.35648326) q[2];
sx q[2];
rz(-2.0429789) q[2];
sx q[2];
rz(-3.0038578) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.288347) q[1];
sx q[1];
rz(-1.4985871) q[1];
sx q[1];
rz(1.9068524) q[1];
rz(0.67948273) q[3];
sx q[3];
rz(-2.3906997) q[3];
sx q[3];
rz(2.7487432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9012458) q[2];
sx q[2];
rz(-2.2525658) q[2];
sx q[2];
rz(-1.3725613) q[2];
rz(1.9339804) q[3];
sx q[3];
rz(-0.6438846) q[3];
sx q[3];
rz(-2.8670368) q[3];
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
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8127301) q[0];
sx q[0];
rz(-2.6565318) q[0];
sx q[0];
rz(3.0364756) q[0];
rz(-2.8016727) q[1];
sx q[1];
rz(-2.2272019) q[1];
sx q[1];
rz(-2.0786659) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.859131) q[0];
sx q[0];
rz(-0.71159186) q[0];
sx q[0];
rz(0.021585394) q[0];
x q[1];
rz(2.379871) q[2];
sx q[2];
rz(-1.3889379) q[2];
sx q[2];
rz(-1.3052502) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.868456) q[1];
sx q[1];
rz(-1.5565762) q[1];
sx q[1];
rz(0.058560024) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.0023554597) q[3];
sx q[3];
rz(-1.071953) q[3];
sx q[3];
rz(2.3209907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.062332705) q[2];
sx q[2];
rz(-1.3881114) q[2];
sx q[2];
rz(1.8049392) q[2];
rz(-3.0873599) q[3];
sx q[3];
rz(-1.1905866) q[3];
sx q[3];
rz(-2.5469053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9688251) q[0];
sx q[0];
rz(-2.4496138) q[0];
sx q[0];
rz(1.2139976) q[0];
rz(-2.0955775) q[1];
sx q[1];
rz(-1.0449301) q[1];
sx q[1];
rz(2.2878343) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19148286) q[0];
sx q[0];
rz(-1.6304509) q[0];
sx q[0];
rz(-1.7050939) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5327024) q[2];
sx q[2];
rz(-0.85431803) q[2];
sx q[2];
rz(-0.68958717) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5484838) q[1];
sx q[1];
rz(-0.9973155) q[1];
sx q[1];
rz(-0.67326633) q[1];
x q[2];
rz(0.3993897) q[3];
sx q[3];
rz(-0.37617427) q[3];
sx q[3];
rz(1.1272421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0115016) q[2];
sx q[2];
rz(-1.6513731) q[2];
sx q[2];
rz(-0.74756527) q[2];
rz(0.93303624) q[3];
sx q[3];
rz(-1.7739762) q[3];
sx q[3];
rz(-0.30195495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1099243) q[0];
sx q[0];
rz(-2.1878991) q[0];
sx q[0];
rz(0.64144301) q[0];
rz(2.172442) q[1];
sx q[1];
rz(-1.0698003) q[1];
sx q[1];
rz(-1.1136805) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4232491) q[0];
sx q[0];
rz(-1.9750496) q[0];
sx q[0];
rz(3.0358311) q[0];
x q[1];
rz(-0.53345726) q[2];
sx q[2];
rz(-1.0048303) q[2];
sx q[2];
rz(0.12803687) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.1195928) q[1];
sx q[1];
rz(-2.5869859) q[1];
sx q[1];
rz(-1.9153992) q[1];
x q[2];
rz(-1.6051172) q[3];
sx q[3];
rz(-2.0908818) q[3];
sx q[3];
rz(-1.8994562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.5668737) q[2];
sx q[2];
rz(-2.4428664) q[2];
sx q[2];
rz(-2.3835772) q[2];
rz(-2.8356683) q[3];
sx q[3];
rz(-0.35861349) q[3];
sx q[3];
rz(0.44736403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4161943) q[0];
sx q[0];
rz(-2.5328126) q[0];
sx q[0];
rz(-0.7533657) q[0];
rz(2.2196409) q[1];
sx q[1];
rz(-1.7148858) q[1];
sx q[1];
rz(3.0070378) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0908703) q[0];
sx q[0];
rz(-0.36733741) q[0];
sx q[0];
rz(-1.814117) q[0];
x q[1];
rz(-2.5002648) q[2];
sx q[2];
rz(-1.3654764) q[2];
sx q[2];
rz(2.5032798) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0363732) q[1];
sx q[1];
rz(-1.4810307) q[1];
sx q[1];
rz(-1.5516993) q[1];
x q[2];
rz(-2.5882072) q[3];
sx q[3];
rz(-2.3883005) q[3];
sx q[3];
rz(-0.91915059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.77384633) q[2];
sx q[2];
rz(-1.6951122) q[2];
sx q[2];
rz(1.9473677) q[2];
rz(-1.9959244) q[3];
sx q[3];
rz(-2.001389) q[3];
sx q[3];
rz(2.0755419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0399748) q[0];
sx q[0];
rz(-2.6823253) q[0];
sx q[0];
rz(2.7591144) q[0];
rz(-0.76599145) q[1];
sx q[1];
rz(-1.6440369) q[1];
sx q[1];
rz(1.8006178) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1857698) q[0];
sx q[0];
rz(-2.9583724) q[0];
sx q[0];
rz(-0.51143666) q[0];
x q[1];
rz(-1.4828311) q[2];
sx q[2];
rz(-0.81552699) q[2];
sx q[2];
rz(1.2620827) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.080440259) q[1];
sx q[1];
rz(-0.78394475) q[1];
sx q[1];
rz(2.817201) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9267335) q[3];
sx q[3];
rz(-1.6481457) q[3];
sx q[3];
rz(-2.429395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0535023) q[2];
sx q[2];
rz(-2.221205) q[2];
sx q[2];
rz(-0.11605334) q[2];
rz(-1.8807489) q[3];
sx q[3];
rz(-2.0033658) q[3];
sx q[3];
rz(-1.6837696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55105227) q[0];
sx q[0];
rz(-3.0278979) q[0];
sx q[0];
rz(-2.9569448) q[0];
rz(2.2231936) q[1];
sx q[1];
rz(-1.5469488) q[1];
sx q[1];
rz(-1.437423) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4686615) q[0];
sx q[0];
rz(-1.7743006) q[0];
sx q[0];
rz(-1.0887906) q[0];
rz(-pi) q[1];
rz(0.080731656) q[2];
sx q[2];
rz(-2.417832) q[2];
sx q[2];
rz(2.2249976) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5823707) q[1];
sx q[1];
rz(-2.5032024) q[1];
sx q[1];
rz(-1.8965782) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4653696) q[3];
sx q[3];
rz(-1.3075324) q[3];
sx q[3];
rz(-0.77419188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.664428) q[2];
sx q[2];
rz(-1.2457341) q[2];
sx q[2];
rz(-2.5028382) q[2];
rz(2.3264558) q[3];
sx q[3];
rz(-2.6705948) q[3];
sx q[3];
rz(0.42069978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4169793) q[0];
sx q[0];
rz(-1.8920349) q[0];
sx q[0];
rz(-2.98988) q[0];
rz(-2.3607415) q[1];
sx q[1];
rz(-1.6897222) q[1];
sx q[1];
rz(2.2471468) q[1];
rz(1.1679756) q[2];
sx q[2];
rz(-0.30277534) q[2];
sx q[2];
rz(1.4390611) q[2];
rz(1.5315957) q[3];
sx q[3];
rz(-1.8146252) q[3];
sx q[3];
rz(-1.6011325) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
