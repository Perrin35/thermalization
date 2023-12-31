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
rz(-0.7601127) q[1];
sx q[1];
rz(-0.7926994) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8108609) q[0];
sx q[0];
rz(-0.90002492) q[0];
sx q[0];
rz(3.1109111) q[0];
rz(-pi) q[1];
rz(2.7723299) q[2];
sx q[2];
rz(-2.5662176) q[2];
sx q[2];
rz(2.0865666) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3789062) q[1];
sx q[1];
rz(-2.0204117) q[1];
sx q[1];
rz(-3.0607037) q[1];
x q[2];
rz(0.80134942) q[3];
sx q[3];
rz(-1.0145463) q[3];
sx q[3];
rz(-0.43503161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4684326) q[2];
sx q[2];
rz(-0.94753733) q[2];
sx q[2];
rz(0.49896487) q[2];
rz(-0.4237825) q[3];
sx q[3];
rz(-1.9215923) q[3];
sx q[3];
rz(3.1288778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.394709) q[0];
sx q[0];
rz(-0.53046662) q[0];
sx q[0];
rz(2.6283426) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5502657) q[2];
sx q[2];
rz(-2.1543243) q[2];
sx q[2];
rz(3.025625) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5181527) q[1];
sx q[1];
rz(-1.232051) q[1];
sx q[1];
rz(-0.63697302) q[1];
rz(-pi) q[2];
rz(-2.9404853) q[3];
sx q[3];
rz(-1.4521986) q[3];
sx q[3];
rz(0.12242854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.83313292) q[2];
sx q[2];
rz(-2.5434727) q[2];
sx q[2];
rz(-1.0773405) q[2];
rz(-0.16263738) q[3];
sx q[3];
rz(-1.2626303) q[3];
sx q[3];
rz(-1.8796273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4804374) q[0];
sx q[0];
rz(-0.85395956) q[0];
sx q[0];
rz(-2.5307632) q[0];
rz(-1.4340596) q[1];
sx q[1];
rz(-1.779498) q[1];
sx q[1];
rz(-2.6909713) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47942802) q[0];
sx q[0];
rz(-0.084708609) q[0];
sx q[0];
rz(-1.1685632) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2497555) q[2];
sx q[2];
rz(-0.46687451) q[2];
sx q[2];
rz(2.4613949) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1994985) q[1];
sx q[1];
rz(-2.4220785) q[1];
sx q[1];
rz(2.662563) q[1];
rz(-pi) q[2];
rz(-1.2080473) q[3];
sx q[3];
rz(-0.72788531) q[3];
sx q[3];
rz(-0.7286275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0314363) q[2];
sx q[2];
rz(-1.0855731) q[2];
sx q[2];
rz(-1.9499367) q[2];
rz(-0.84447652) q[3];
sx q[3];
rz(-2.8561487) q[3];
sx q[3];
rz(-0.58117956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0464756) q[0];
sx q[0];
rz(-1.7549544) q[0];
sx q[0];
rz(-1.1412096) q[0];
rz(1.8449239) q[1];
sx q[1];
rz(-2.7099687) q[1];
sx q[1];
rz(2.8362714) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5019708) q[0];
sx q[0];
rz(-1.0143711) q[0];
sx q[0];
rz(1.200838) q[0];
x q[1];
rz(1.9919473) q[2];
sx q[2];
rz(-1.9990168) q[2];
sx q[2];
rz(-0.8840094) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.893559) q[1];
sx q[1];
rz(-2.3602924) q[1];
sx q[1];
rz(-1.7209189) q[1];
x q[2];
rz(-1.0563072) q[3];
sx q[3];
rz(-1.1122276) q[3];
sx q[3];
rz(1.4460627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9332463) q[2];
sx q[2];
rz(-1.5758005) q[2];
sx q[2];
rz(2.7944881) q[2];
rz(2.7582205) q[3];
sx q[3];
rz(-0.88630668) q[3];
sx q[3];
rz(1.0679519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7326996) q[0];
sx q[0];
rz(-2.2786031) q[0];
sx q[0];
rz(0.53891671) q[0];
rz(1.7319038) q[1];
sx q[1];
rz(-2.6795645) q[1];
sx q[1];
rz(2.7358823) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6255394) q[0];
sx q[0];
rz(-2.2279983) q[0];
sx q[0];
rz(-1.2217193) q[0];
x q[1];
rz(1.97498) q[2];
sx q[2];
rz(-0.6037854) q[2];
sx q[2];
rz(-3.0886138) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.09771422) q[1];
sx q[1];
rz(-1.8227302) q[1];
sx q[1];
rz(-1.9240727) q[1];
rz(-1.3375072) q[3];
sx q[3];
rz(-0.81987112) q[3];
sx q[3];
rz(1.6398721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.47200176) q[2];
sx q[2];
rz(-2.3386023) q[2];
sx q[2];
rz(-1.6873138) q[2];
rz(0.83868319) q[3];
sx q[3];
rz(-1.2861774) q[3];
sx q[3];
rz(1.3735166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31081653) q[0];
sx q[0];
rz(-0.61284471) q[0];
sx q[0];
rz(-0.40711656) q[0];
rz(0.85917568) q[1];
sx q[1];
rz(-1.5770864) q[1];
sx q[1];
rz(-1.3173332) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2989527) q[0];
sx q[0];
rz(-1.5931804) q[0];
sx q[0];
rz(-1.0277261) q[0];
rz(-pi) q[1];
rz(0.69465722) q[2];
sx q[2];
rz(-0.79840556) q[2];
sx q[2];
rz(-0.77381221) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0560318) q[1];
sx q[1];
rz(-2.3059658) q[1];
sx q[1];
rz(3.1378156) q[1];
rz(0.96885724) q[3];
sx q[3];
rz(-0.67776206) q[3];
sx q[3];
rz(-0.70167752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.09981) q[2];
sx q[2];
rz(-0.87974292) q[2];
sx q[2];
rz(2.2463634) q[2];
rz(-1.2081395) q[3];
sx q[3];
rz(-2.4619305) q[3];
sx q[3];
rz(-0.7081379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5773425) q[0];
sx q[0];
rz(-2.0944493) q[0];
sx q[0];
rz(2.0136925) q[0];
rz(-0.212542) q[1];
sx q[1];
rz(-1.3393341) q[1];
sx q[1];
rz(-1.9099265) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8679778) q[0];
sx q[0];
rz(-0.6379188) q[0];
sx q[0];
rz(0.62423737) q[0];
rz(0.64078757) q[2];
sx q[2];
rz(-2.1038342) q[2];
sx q[2];
rz(1.8719045) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.54214189) q[1];
sx q[1];
rz(-0.25240024) q[1];
sx q[1];
rz(-2.7155994) q[1];
rz(-pi) q[2];
rz(2.1752862) q[3];
sx q[3];
rz(-1.55282) q[3];
sx q[3];
rz(3.0997961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.82905519) q[2];
sx q[2];
rz(-2.2237491) q[2];
sx q[2];
rz(0.60116872) q[2];
rz(0.51856315) q[3];
sx q[3];
rz(-0.31373599) q[3];
sx q[3];
rz(1.4275838) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9940051) q[0];
sx q[0];
rz(-0.00031358263) q[0];
sx q[0];
rz(-3.0684379) q[0];
rz(1.0271429) q[1];
sx q[1];
rz(-1.606753) q[1];
sx q[1];
rz(-0.59246078) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6153529) q[0];
sx q[0];
rz(-0.59229367) q[0];
sx q[0];
rz(-3.0583529) q[0];
rz(-pi) q[1];
x q[1];
rz(0.57581298) q[2];
sx q[2];
rz(-1.1278707) q[2];
sx q[2];
rz(-1.5695614) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2295099) q[1];
sx q[1];
rz(-2.3778937) q[1];
sx q[1];
rz(-2.0565226) q[1];
x q[2];
rz(1.1660277) q[3];
sx q[3];
rz(-1.3169857) q[3];
sx q[3];
rz(0.39241957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.56117326) q[2];
sx q[2];
rz(-1.6757123) q[2];
sx q[2];
rz(-0.44211659) q[2];
rz(-0.90028611) q[3];
sx q[3];
rz(-2.0787652) q[3];
sx q[3];
rz(-3.1282848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(0.028037926) q[0];
sx q[0];
rz(-2.2224764) q[0];
sx q[0];
rz(-1.8369209) q[0];
rz(-1.4470944) q[1];
sx q[1];
rz(-1.1728975) q[1];
sx q[1];
rz(-0.56217271) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2660414) q[0];
sx q[0];
rz(-2.3444359) q[0];
sx q[0];
rz(-2.6334727) q[0];
x q[1];
rz(-2.1015342) q[2];
sx q[2];
rz(-1.6518403) q[2];
sx q[2];
rz(-2.9329252) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.11645424) q[1];
sx q[1];
rz(-1.0944195) q[1];
sx q[1];
rz(-2.3548747) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.4060679) q[3];
sx q[3];
rz(-1.1998402) q[3];
sx q[3];
rz(-3.1379116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.28329453) q[2];
sx q[2];
rz(-1.4026182) q[2];
sx q[2];
rz(-0.98158681) q[2];
rz(-0.0043407241) q[3];
sx q[3];
rz(-2.0948295) q[3];
sx q[3];
rz(0.59345746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8269862) q[0];
sx q[0];
rz(-1.9208603) q[0];
sx q[0];
rz(0.04709588) q[0];
rz(1.3866407) q[1];
sx q[1];
rz(-1.1562693) q[1];
sx q[1];
rz(0.047853619) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3536516) q[0];
sx q[0];
rz(-1.5418058) q[0];
sx q[0];
rz(2.2257462) q[0];
rz(-pi) q[1];
rz(3.0907862) q[2];
sx q[2];
rz(-2.0127957) q[2];
sx q[2];
rz(2.6450784) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0117482) q[1];
sx q[1];
rz(-2.1455892) q[1];
sx q[1];
rz(0.5188491) q[1];
rz(-1.5157394) q[3];
sx q[3];
rz(-1.1686472) q[3];
sx q[3];
rz(-2.1189342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.64426595) q[2];
sx q[2];
rz(-1.9571807) q[2];
sx q[2];
rz(-2.2052374) q[2];
rz(0.80091536) q[3];
sx q[3];
rz(-2.622486) q[3];
sx q[3];
rz(2.463533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3180278) q[0];
sx q[0];
rz(-2.2021273) q[0];
sx q[0];
rz(0.20903023) q[0];
rz(-2.6387852) q[1];
sx q[1];
rz(-1.8223169) q[1];
sx q[1];
rz(1.0261789) q[1];
rz(-0.55124333) q[2];
sx q[2];
rz(-1.5521282) q[2];
sx q[2];
rz(2.9350401) q[2];
rz(-2.3021163) q[3];
sx q[3];
rz(-0.52121938) q[3];
sx q[3];
rz(-2.2073707) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
