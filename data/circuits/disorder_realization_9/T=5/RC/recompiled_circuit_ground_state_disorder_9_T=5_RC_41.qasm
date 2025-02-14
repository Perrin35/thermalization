OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.61152148) q[0];
sx q[0];
rz(5.2103142) q[0];
sx q[0];
rz(11.13902) q[0];
rz(0.99675769) q[1];
sx q[1];
rz(-2.5280894) q[1];
sx q[1];
rz(0.065453425) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6320541) q[0];
sx q[0];
rz(-0.51842123) q[0];
sx q[0];
rz(-2.753756) q[0];
rz(-pi) q[1];
x q[1];
rz(0.10721389) q[2];
sx q[2];
rz(-0.37069025) q[2];
sx q[2];
rz(3.0840741) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.10852818) q[1];
sx q[1];
rz(-2.4208596) q[1];
sx q[1];
rz(-1.1041376) q[1];
rz(2.5259326) q[3];
sx q[3];
rz(-1.1634852) q[3];
sx q[3];
rz(-0.034192702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2685711) q[2];
sx q[2];
rz(-2.7950588) q[2];
sx q[2];
rz(1.2968501) q[2];
rz(-1.1067363) q[3];
sx q[3];
rz(-0.96702558) q[3];
sx q[3];
rz(1.6905674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1065555) q[0];
sx q[0];
rz(-0.7283926) q[0];
sx q[0];
rz(-1.534071) q[0];
rz(-2.2976177) q[1];
sx q[1];
rz(-1.5183828) q[1];
sx q[1];
rz(-0.27110505) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8753238) q[0];
sx q[0];
rz(-1.3014587) q[0];
sx q[0];
rz(1.6104417) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9366802) q[2];
sx q[2];
rz(-0.26981631) q[2];
sx q[2];
rz(-1.7797949) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0738759) q[1];
sx q[1];
rz(-0.42781204) q[1];
sx q[1];
rz(-2.907987) q[1];
rz(1.4976147) q[3];
sx q[3];
rz(-1.7203619) q[3];
sx q[3];
rz(2.4084546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8035182) q[2];
sx q[2];
rz(-2.1536571) q[2];
sx q[2];
rz(0.79338497) q[2];
rz(3.0986943) q[3];
sx q[3];
rz(-1.2267313) q[3];
sx q[3];
rz(2.4734917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-1.0637829) q[0];
sx q[0];
rz(-0.48352799) q[0];
sx q[0];
rz(-3.0384645) q[0];
rz(-2.8887796) q[1];
sx q[1];
rz(-1.4272855) q[1];
sx q[1];
rz(2.8275183) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47427705) q[0];
sx q[0];
rz(-1.0768824) q[0];
sx q[0];
rz(2.0163573) q[0];
rz(-pi) q[1];
rz(0.057217802) q[2];
sx q[2];
rz(-2.0524745) q[2];
sx q[2];
rz(1.526598) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9772226) q[1];
sx q[1];
rz(-0.80406351) q[1];
sx q[1];
rz(-2.3146446) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5669085) q[3];
sx q[3];
rz(-2.8369224) q[3];
sx q[3];
rz(-1.6449354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.126943) q[2];
sx q[2];
rz(-2.5003771) q[2];
sx q[2];
rz(-2.3250697) q[2];
rz(0.71216312) q[3];
sx q[3];
rz(-1.687259) q[3];
sx q[3];
rz(-0.85675353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6461058) q[0];
sx q[0];
rz(-3.0480338) q[0];
sx q[0];
rz(-0.58864546) q[0];
rz(-0.37216392) q[1];
sx q[1];
rz(-1.2969505) q[1];
sx q[1];
rz(-0.94474307) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0030768) q[0];
sx q[0];
rz(-2.8154065) q[0];
sx q[0];
rz(-1.7762) q[0];
rz(-pi) q[1];
rz(2.0633374) q[2];
sx q[2];
rz(-1.8027935) q[2];
sx q[2];
rz(1.3859232) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1334298) q[1];
sx q[1];
rz(-1.9979229) q[1];
sx q[1];
rz(-2.5168946) q[1];
rz(1.9362336) q[3];
sx q[3];
rz(-1.7492709) q[3];
sx q[3];
rz(0.5040666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9113691) q[2];
sx q[2];
rz(-1.5350124) q[2];
sx q[2];
rz(1.4851419) q[2];
rz(-0.39561513) q[3];
sx q[3];
rz(-1.4916568) q[3];
sx q[3];
rz(3.0806471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38451251) q[0];
sx q[0];
rz(-2.7482996) q[0];
sx q[0];
rz(1.242189) q[0];
rz(-0.077266129) q[1];
sx q[1];
rz(-1.9237513) q[1];
sx q[1];
rz(-0.088931106) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9327527) q[0];
sx q[0];
rz(-1.7249301) q[0];
sx q[0];
rz(-0.16618118) q[0];
rz(-0.16392608) q[2];
sx q[2];
rz(-0.55183119) q[2];
sx q[2];
rz(1.1140119) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4684852) q[1];
sx q[1];
rz(-1.2947646) q[1];
sx q[1];
rz(-2.0620146) q[1];
rz(2.2054178) q[3];
sx q[3];
rz(-0.16491865) q[3];
sx q[3];
rz(2.0009964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.58101216) q[2];
sx q[2];
rz(-1.7580914) q[2];
sx q[2];
rz(-1.9182473) q[2];
rz(1.9129725) q[3];
sx q[3];
rz(-0.88594985) q[3];
sx q[3];
rz(-2.3942053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0864047) q[0];
sx q[0];
rz(-2.3265657) q[0];
sx q[0];
rz(-2.1888457) q[0];
rz(-1.791026) q[1];
sx q[1];
rz(-1.2341276) q[1];
sx q[1];
rz(-1.9780673) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.945283) q[0];
sx q[0];
rz(-2.3988814) q[0];
sx q[0];
rz(0.16603827) q[0];
x q[1];
rz(2.7609652) q[2];
sx q[2];
rz(-0.63063313) q[2];
sx q[2];
rz(-2.6826114) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.42683329) q[1];
sx q[1];
rz(-1.3275255) q[1];
sx q[1];
rz(-0.053364874) q[1];
rz(-pi) q[2];
rz(2.5599078) q[3];
sx q[3];
rz(-0.69184985) q[3];
sx q[3];
rz(-2.8210616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.20840883) q[2];
sx q[2];
rz(-1.2196502) q[2];
sx q[2];
rz(-1.5677933) q[2];
rz(-0.21582223) q[3];
sx q[3];
rz(-0.62041557) q[3];
sx q[3];
rz(-0.1568493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.02483524) q[0];
sx q[0];
rz(-0.76501608) q[0];
sx q[0];
rz(-1.5417954) q[0];
rz(2.722591) q[1];
sx q[1];
rz(-1.1083138) q[1];
sx q[1];
rz(1.3655183) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5758662) q[0];
sx q[0];
rz(-0.64896261) q[0];
sx q[0];
rz(-1.2470232) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.55115564) q[2];
sx q[2];
rz(-0.86116448) q[2];
sx q[2];
rz(0.9786419) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8518062) q[1];
sx q[1];
rz(-0.86871457) q[1];
sx q[1];
rz(-2.9243824) q[1];
x q[2];
rz(-2.0715874) q[3];
sx q[3];
rz(-2.4907101) q[3];
sx q[3];
rz(0.18903519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2447723) q[2];
sx q[2];
rz(-3.012487) q[2];
sx q[2];
rz(1.6884035) q[2];
rz(-1.7425273) q[3];
sx q[3];
rz(-1.1996256) q[3];
sx q[3];
rz(2.3340268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6145265) q[0];
sx q[0];
rz(-2.2699321) q[0];
sx q[0];
rz(2.6283619) q[0];
rz(-2.1649583) q[1];
sx q[1];
rz(-2.4726424) q[1];
sx q[1];
rz(1.4248779) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7920765) q[0];
sx q[0];
rz(-0.87679014) q[0];
sx q[0];
rz(-0.16887878) q[0];
rz(-0.54855697) q[2];
sx q[2];
rz(-1.6335677) q[2];
sx q[2];
rz(1.5323764) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1816563) q[1];
sx q[1];
rz(-1.8605369) q[1];
sx q[1];
rz(-1.1162037) q[1];
rz(-pi) q[2];
rz(0.11792169) q[3];
sx q[3];
rz(-1.7365082) q[3];
sx q[3];
rz(-2.5664701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6500515) q[2];
sx q[2];
rz(-1.0074793) q[2];
sx q[2];
rz(0.36637351) q[2];
rz(2.5036687) q[3];
sx q[3];
rz(-1.7584691) q[3];
sx q[3];
rz(-1.5095507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.76846182) q[0];
sx q[0];
rz(-2.8167384) q[0];
sx q[0];
rz(1.984206) q[0];
rz(2.6034082) q[1];
sx q[1];
rz(-1.3342074) q[1];
sx q[1];
rz(-3.11788) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2185604) q[0];
sx q[0];
rz(-2.5171748) q[0];
sx q[0];
rz(0.10186457) q[0];
rz(-pi) q[1];
rz(-1.1309409) q[2];
sx q[2];
rz(-0.80200101) q[2];
sx q[2];
rz(1.8126876) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4166322) q[1];
sx q[1];
rz(-1.2051996) q[1];
sx q[1];
rz(-2.599692) q[1];
rz(-pi) q[2];
rz(1.7718763) q[3];
sx q[3];
rz(-1.4477125) q[3];
sx q[3];
rz(-1.1986002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9862765) q[2];
sx q[2];
rz(-1.3331058) q[2];
sx q[2];
rz(-2.7434366) q[2];
rz(0.46477535) q[3];
sx q[3];
rz(-2.0538752) q[3];
sx q[3];
rz(1.3409415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.583113) q[0];
sx q[0];
rz(-0.3322424) q[0];
sx q[0];
rz(0.94859052) q[0];
rz(0.96725431) q[1];
sx q[1];
rz(-1.9890669) q[1];
sx q[1];
rz(-2.1327877) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0111158) q[0];
sx q[0];
rz(-1.563894) q[0];
sx q[0];
rz(-2.3409977) q[0];
rz(-2.45935) q[2];
sx q[2];
rz(-0.80305225) q[2];
sx q[2];
rz(-2.19953) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.12647835) q[1];
sx q[1];
rz(-2.6228944) q[1];
sx q[1];
rz(2.7706503) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7901374) q[3];
sx q[3];
rz(-0.92445856) q[3];
sx q[3];
rz(-1.5953894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.41445109) q[2];
sx q[2];
rz(-2.264617) q[2];
sx q[2];
rz(0.021075185) q[2];
rz(0.4840788) q[3];
sx q[3];
rz(-2.0120072) q[3];
sx q[3];
rz(2.4626203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2912343) q[0];
sx q[0];
rz(-1.6073011) q[0];
sx q[0];
rz(-1.4144443) q[0];
rz(2.5454632) q[1];
sx q[1];
rz(-0.78501751) q[1];
sx q[1];
rz(2.6529978) q[1];
rz(0.1487371) q[2];
sx q[2];
rz(-1.7784468) q[2];
sx q[2];
rz(-0.39685984) q[2];
rz(-1.7763185) q[3];
sx q[3];
rz(-1.4649656) q[3];
sx q[3];
rz(-3.0595915) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
