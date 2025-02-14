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
rz(-0.64033163) q[0];
sx q[0];
rz(-0.92136541) q[0];
sx q[0];
rz(1.6093572) q[0];
rz(0.20707239) q[1];
sx q[1];
rz(4.3011811) q[1];
sx q[1];
rz(11.094697) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64016104) q[0];
sx q[0];
rz(-0.45217538) q[0];
sx q[0];
rz(-0.78340952) q[0];
x q[1];
rz(0.70090909) q[2];
sx q[2];
rz(-0.94628382) q[2];
sx q[2];
rz(-2.8025742) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0826559) q[1];
sx q[1];
rz(-1.3322222) q[1];
sx q[1];
rz(0.10175565) q[1];
rz(-1.3892646) q[3];
sx q[3];
rz(-2.1588491) q[3];
sx q[3];
rz(-1.4168038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6472935) q[2];
sx q[2];
rz(-2.1493201) q[2];
sx q[2];
rz(2.5956019) q[2];
rz(-0.93849385) q[3];
sx q[3];
rz(-2.6280845) q[3];
sx q[3];
rz(0.45309666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12478011) q[0];
sx q[0];
rz(-1.1174959) q[0];
sx q[0];
rz(0.6413396) q[0];
rz(-1.9524139) q[1];
sx q[1];
rz(-0.42354217) q[1];
sx q[1];
rz(-1.10434) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1416407) q[0];
sx q[0];
rz(-1.4237836) q[0];
sx q[0];
rz(2.5070058) q[0];
rz(2.10675) q[2];
sx q[2];
rz(-1.8551197) q[2];
sx q[2];
rz(-0.51229561) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.84766372) q[1];
sx q[1];
rz(-1.0231471) q[1];
sx q[1];
rz(2.6002162) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5303753) q[3];
sx q[3];
rz(-1.6580402) q[3];
sx q[3];
rz(-0.62007346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.39893338) q[2];
sx q[2];
rz(-2.0519966) q[2];
sx q[2];
rz(0.99417865) q[2];
rz(1.582877) q[3];
sx q[3];
rz(-2.4623058) q[3];
sx q[3];
rz(-2.5905632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-3.1347374) q[0];
sx q[0];
rz(-2.0515433) q[0];
sx q[0];
rz(2.5110974) q[0];
rz(-0.14671239) q[1];
sx q[1];
rz(-1.881733) q[1];
sx q[1];
rz(2.2781118) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9778507) q[0];
sx q[0];
rz(-1.2702748) q[0];
sx q[0];
rz(-2.252821) q[0];
x q[1];
rz(-1.2892154) q[2];
sx q[2];
rz(-1.3669434) q[2];
sx q[2];
rz(1.8903763) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.32541986) q[1];
sx q[1];
rz(-2.6920941) q[1];
sx q[1];
rz(-1.3260496) q[1];
rz(-pi) q[2];
rz(-2.9894838) q[3];
sx q[3];
rz(-0.31732355) q[3];
sx q[3];
rz(-0.0025686669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8063987) q[2];
sx q[2];
rz(-1.4525745) q[2];
sx q[2];
rz(0.88649583) q[2];
rz(2.7116306) q[3];
sx q[3];
rz(-2.6468247) q[3];
sx q[3];
rz(-0.89216843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48915136) q[0];
sx q[0];
rz(-0.2944856) q[0];
sx q[0];
rz(0.44892204) q[0];
rz(2.4849675) q[1];
sx q[1];
rz(-1.4337599) q[1];
sx q[1];
rz(2.8112559) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93189222) q[0];
sx q[0];
rz(-2.7468514) q[0];
sx q[0];
rz(2.1861211) q[0];
x q[1];
rz(-2.0520211) q[2];
sx q[2];
rz(-0.46483609) q[2];
sx q[2];
rz(-0.46457738) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9724839) q[1];
sx q[1];
rz(-2.0044998) q[1];
sx q[1];
rz(0.32167158) q[1];
rz(-1.6918332) q[3];
sx q[3];
rz(-1.5618298) q[3];
sx q[3];
rz(-3.1154273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0535447) q[2];
sx q[2];
rz(-1.6036754) q[2];
sx q[2];
rz(0.21928445) q[2];
rz(2.3413279) q[3];
sx q[3];
rz(-2.181668) q[3];
sx q[3];
rz(0.76103359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6680172) q[0];
sx q[0];
rz(-2.0125084) q[0];
sx q[0];
rz(-1.6582723) q[0];
rz(-2.5738916) q[1];
sx q[1];
rz(-1.2115819) q[1];
sx q[1];
rz(1.6426881) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9358112) q[0];
sx q[0];
rz(-1.0438127) q[0];
sx q[0];
rz(-1.7409659) q[0];
rz(-pi) q[1];
rz(-0.11168555) q[2];
sx q[2];
rz(-3.0431261) q[2];
sx q[2];
rz(0.51375997) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.066047) q[1];
sx q[1];
rz(-2.4986364) q[1];
sx q[1];
rz(-2.8075904) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.12171774) q[3];
sx q[3];
rz(-1.7608402) q[3];
sx q[3];
rz(-2.7167883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6296926) q[2];
sx q[2];
rz(-2.3834159) q[2];
sx q[2];
rz(-2.8153815) q[2];
rz(-0.1263667) q[3];
sx q[3];
rz(-2.5765403) q[3];
sx q[3];
rz(2.8700184) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3417974) q[0];
sx q[0];
rz(-1.3038776) q[0];
sx q[0];
rz(2.6190992) q[0];
rz(-2.3937468) q[1];
sx q[1];
rz(-1.5380842) q[1];
sx q[1];
rz(1.1572256) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8337473) q[0];
sx q[0];
rz(-0.54824775) q[0];
sx q[0];
rz(-0.98497106) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0339478) q[2];
sx q[2];
rz(-1.2105296) q[2];
sx q[2];
rz(2.2778794) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6023691) q[1];
sx q[1];
rz(-1.4280978) q[1];
sx q[1];
rz(-3.0485832) q[1];
rz(-pi) q[2];
x q[2];
rz(0.32706599) q[3];
sx q[3];
rz(-1.4983791) q[3];
sx q[3];
rz(-1.678997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.130827) q[2];
sx q[2];
rz(-0.5126493) q[2];
sx q[2];
rz(1.1831247) q[2];
rz(3.0986687) q[3];
sx q[3];
rz(-1.639651) q[3];
sx q[3];
rz(-0.24421282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44508988) q[0];
sx q[0];
rz(-2.255991) q[0];
sx q[0];
rz(-0.67286277) q[0];
rz(-2.591835) q[1];
sx q[1];
rz(-1.8472698) q[1];
sx q[1];
rz(-1.4901644) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5527809) q[0];
sx q[0];
rz(-2.5500265) q[0];
sx q[0];
rz(0.32365303) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3734858) q[2];
sx q[2];
rz(-0.66181493) q[2];
sx q[2];
rz(2.2827374) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.71201176) q[1];
sx q[1];
rz(-2.6954571) q[1];
sx q[1];
rz(2.9175116) q[1];
rz(-2.1792322) q[3];
sx q[3];
rz(-0.058789805) q[3];
sx q[3];
rz(0.99942452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0979536) q[2];
sx q[2];
rz(-2.592228) q[2];
sx q[2];
rz(-1.1489541) q[2];
rz(-2.7502381) q[3];
sx q[3];
rz(-1.9214182) q[3];
sx q[3];
rz(-2.2933188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-1.5452165) q[0];
sx q[0];
rz(-1.9356118) q[0];
sx q[0];
rz(0.69860953) q[0];
rz(-0.85083234) q[1];
sx q[1];
rz(-2.5406676) q[1];
sx q[1];
rz(-0.94815475) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8671001) q[0];
sx q[0];
rz(-0.72625181) q[0];
sx q[0];
rz(-0.12402835) q[0];
x q[1];
rz(0.053669503) q[2];
sx q[2];
rz(-1.9000152) q[2];
sx q[2];
rz(-1.7754796) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8906456) q[1];
sx q[1];
rz(-1.0140103) q[1];
sx q[1];
rz(0.48018881) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.83017577) q[3];
sx q[3];
rz(-2.1599033) q[3];
sx q[3];
rz(2.0996706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4725388) q[2];
sx q[2];
rz(-1.6927745) q[2];
sx q[2];
rz(-1.4581468) q[2];
rz(-2.0504047) q[3];
sx q[3];
rz(-0.89717054) q[3];
sx q[3];
rz(0.18479656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75337306) q[0];
sx q[0];
rz(-0.18167697) q[0];
sx q[0];
rz(1.3004119) q[0];
rz(0.45937195) q[1];
sx q[1];
rz(-1.4540902) q[1];
sx q[1];
rz(-0.5836817) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4874909) q[0];
sx q[0];
rz(-1.0922474) q[0];
sx q[0];
rz(0.74639456) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.74588512) q[2];
sx q[2];
rz(-3.0891232) q[2];
sx q[2];
rz(1.3191163) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.65185279) q[1];
sx q[1];
rz(-1.0612584) q[1];
sx q[1];
rz(-2.7368828) q[1];
rz(-pi) q[2];
rz(-1.0631869) q[3];
sx q[3];
rz(-0.48590966) q[3];
sx q[3];
rz(-2.8775281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5251081) q[2];
sx q[2];
rz(-2.8664092) q[2];
sx q[2];
rz(-2.6194561) q[2];
rz(2.3790242) q[3];
sx q[3];
rz(-1.6430166) q[3];
sx q[3];
rz(0.34408072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6793215) q[0];
sx q[0];
rz(-1.2079879) q[0];
sx q[0];
rz(-0.59360498) q[0];
rz(2.4484334) q[1];
sx q[1];
rz(-0.79265541) q[1];
sx q[1];
rz(-0.75137442) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5042413) q[0];
sx q[0];
rz(-2.6587186) q[0];
sx q[0];
rz(0.58321799) q[0];
rz(0.85752731) q[2];
sx q[2];
rz(-0.81999841) q[2];
sx q[2];
rz(-1.6423026) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.33022949) q[1];
sx q[1];
rz(-2.4733481) q[1];
sx q[1];
rz(-0.10295566) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4937711) q[3];
sx q[3];
rz(-1.739177) q[3];
sx q[3];
rz(-1.8986957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1410602) q[2];
sx q[2];
rz(-2.4324721) q[2];
sx q[2];
rz(3.1179023) q[2];
rz(1.116811) q[3];
sx q[3];
rz(-1.4477372) q[3];
sx q[3];
rz(2.007133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5524207) q[0];
sx q[0];
rz(-1.9235274) q[0];
sx q[0];
rz(-2.0148475) q[0];
rz(1.8213656) q[1];
sx q[1];
rz(-1.2757433) q[1];
sx q[1];
rz(-2.129) q[1];
rz(-2.597328) q[2];
sx q[2];
rz(-1.7901292) q[2];
sx q[2];
rz(-2.3397056) q[2];
rz(-0.79441073) q[3];
sx q[3];
rz(-2.8230739) q[3];
sx q[3];
rz(2.4369797) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
