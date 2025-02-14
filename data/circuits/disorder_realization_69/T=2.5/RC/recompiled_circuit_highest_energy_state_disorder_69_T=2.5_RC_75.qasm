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
rz(1.5975098) q[0];
sx q[0];
rz(-1.37356) q[0];
sx q[0];
rz(-2.1231667) q[0];
rz(2.6234558) q[1];
sx q[1];
rz(-0.75704804) q[1];
sx q[1];
rz(2.5098324) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1904097) q[0];
sx q[0];
rz(-2.8993239) q[0];
sx q[0];
rz(2.1184068) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4912946) q[2];
sx q[2];
rz(-0.30361807) q[2];
sx q[2];
rz(-2.7111766) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.5093928) q[1];
sx q[1];
rz(-2.191847) q[1];
sx q[1];
rz(-0.34624798) q[1];
rz(-pi) q[2];
rz(-1.8970892) q[3];
sx q[3];
rz(-1.658381) q[3];
sx q[3];
rz(1.0911986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.53039256) q[2];
sx q[2];
rz(-2.7860614) q[2];
sx q[2];
rz(1.4872888) q[2];
rz(-0.90855956) q[3];
sx q[3];
rz(-0.23659758) q[3];
sx q[3];
rz(1.0049741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0205883) q[0];
sx q[0];
rz(-1.4326743) q[0];
sx q[0];
rz(2.9793136) q[0];
rz(-0.24457112) q[1];
sx q[1];
rz(-1.2030315) q[1];
sx q[1];
rz(0.29168209) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26839089) q[0];
sx q[0];
rz(-1.8531728) q[0];
sx q[0];
rz(-1.3703521) q[0];
rz(2.2333916) q[2];
sx q[2];
rz(-0.35396265) q[2];
sx q[2];
rz(-0.56828729) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.32714601) q[1];
sx q[1];
rz(-1.8015878) q[1];
sx q[1];
rz(0.77951851) q[1];
rz(-pi) q[2];
rz(-0.41145153) q[3];
sx q[3];
rz(-1.9570159) q[3];
sx q[3];
rz(-0.4377818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.67819277) q[2];
sx q[2];
rz(-0.86898154) q[2];
sx q[2];
rz(-1.0758859) q[2];
rz(1.2470657) q[3];
sx q[3];
rz(-1.5970634) q[3];
sx q[3];
rz(-1.4472848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70188824) q[0];
sx q[0];
rz(-1.1264369) q[0];
sx q[0];
rz(0.18371789) q[0];
rz(1.5048997) q[1];
sx q[1];
rz(-2.3972062) q[1];
sx q[1];
rz(2.9768129) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4713584) q[0];
sx q[0];
rz(-0.58588282) q[0];
sx q[0];
rz(1.5878994) q[0];
rz(-pi) q[1];
x q[1];
rz(0.38305958) q[2];
sx q[2];
rz(-2.428741) q[2];
sx q[2];
rz(1.0733611) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2307869) q[1];
sx q[1];
rz(-1.1890829) q[1];
sx q[1];
rz(0.64818212) q[1];
x q[2];
rz(-2.80071) q[3];
sx q[3];
rz(-1.0629144) q[3];
sx q[3];
rz(1.9991877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.034721) q[2];
sx q[2];
rz(-1.160459) q[2];
sx q[2];
rz(1.1064233) q[2];
rz(-0.70118457) q[3];
sx q[3];
rz(-1.3283575) q[3];
sx q[3];
rz(2.6760694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3727386) q[0];
sx q[0];
rz(-1.1356069) q[0];
sx q[0];
rz(0.94451529) q[0];
rz(-1.5185897) q[1];
sx q[1];
rz(-0.98097643) q[1];
sx q[1];
rz(1.6015582) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.042451579) q[0];
sx q[0];
rz(-0.97223982) q[0];
sx q[0];
rz(-2.4392946) q[0];
rz(-1.6883019) q[2];
sx q[2];
rz(-1.4521993) q[2];
sx q[2];
rz(2.8343253) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.67805144) q[1];
sx q[1];
rz(-1.1992595) q[1];
sx q[1];
rz(-2.3080565) q[1];
rz(-1.6155195) q[3];
sx q[3];
rz(-0.8208684) q[3];
sx q[3];
rz(-0.57425153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.010178415) q[2];
sx q[2];
rz(-2.5203036) q[2];
sx q[2];
rz(-2.9739001) q[2];
rz(-0.0532648) q[3];
sx q[3];
rz(-1.1054509) q[3];
sx q[3];
rz(1.7648599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.57529706) q[0];
sx q[0];
rz(-0.10558858) q[0];
sx q[0];
rz(0.29933023) q[0];
rz(-2.7686367) q[1];
sx q[1];
rz(-1.8920218) q[1];
sx q[1];
rz(1.6711055) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1363935) q[0];
sx q[0];
rz(-1.0233425) q[0];
sx q[0];
rz(1.0639079) q[0];
rz(-pi) q[1];
rz(-2.6777381) q[2];
sx q[2];
rz(-1.6783444) q[2];
sx q[2];
rz(0.35075089) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8052696) q[1];
sx q[1];
rz(-1.0674369) q[1];
sx q[1];
rz(-2.4783647) q[1];
rz(-pi) q[2];
rz(2.1924952) q[3];
sx q[3];
rz(-1.7142363) q[3];
sx q[3];
rz(-0.52191097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2028929) q[2];
sx q[2];
rz(-0.6627658) q[2];
sx q[2];
rz(1.1104442) q[2];
rz(-2.2796196) q[3];
sx q[3];
rz(-1.5098666) q[3];
sx q[3];
rz(1.8814258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92457572) q[0];
sx q[0];
rz(-0.68704263) q[0];
sx q[0];
rz(-2.7365015) q[0];
rz(0.336126) q[1];
sx q[1];
rz(-1.4210217) q[1];
sx q[1];
rz(-2.1902671) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38430957) q[0];
sx q[0];
rz(-1.8258137) q[0];
sx q[0];
rz(2.7268305) q[0];
rz(-1.3589922) q[2];
sx q[2];
rz(-0.84583218) q[2];
sx q[2];
rz(2.447809) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3094294) q[1];
sx q[1];
rz(-2.1346666) q[1];
sx q[1];
rz(1.840074) q[1];
rz(-2.5441465) q[3];
sx q[3];
rz(-2.8077112) q[3];
sx q[3];
rz(-2.2528439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.090791) q[2];
sx q[2];
rz(-1.7369221) q[2];
sx q[2];
rz(0.23078272) q[2];
rz(2.9595621) q[3];
sx q[3];
rz(-2.5735276) q[3];
sx q[3];
rz(-0.52869421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7886605) q[0];
sx q[0];
rz(-0.039529888) q[0];
sx q[0];
rz(-0.93210644) q[0];
rz(0.63356361) q[1];
sx q[1];
rz(-2.0220058) q[1];
sx q[1];
rz(-1.8036141) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.103823) q[0];
sx q[0];
rz(-1.4892206) q[0];
sx q[0];
rz(-1.6768811) q[0];
rz(-pi) q[1];
x q[1];
rz(0.39057486) q[2];
sx q[2];
rz(-2.8468067) q[2];
sx q[2];
rz(0.24402555) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7576068) q[1];
sx q[1];
rz(-1.9644871) q[1];
sx q[1];
rz(2.2298953) q[1];
x q[2];
rz(0.49832817) q[3];
sx q[3];
rz(-1.5428203) q[3];
sx q[3];
rz(1.2714112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6771217) q[2];
sx q[2];
rz(-1.9147583) q[2];
sx q[2];
rz(-1.202549) q[2];
rz(-0.36457148) q[3];
sx q[3];
rz(-0.69260827) q[3];
sx q[3];
rz(2.306849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4778336) q[0];
sx q[0];
rz(-0.57357016) q[0];
sx q[0];
rz(-2.9649576) q[0];
rz(-0.44003507) q[1];
sx q[1];
rz(-1.1853848) q[1];
sx q[1];
rz(0.96493351) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4266708) q[0];
sx q[0];
rz(-1.7650801) q[0];
sx q[0];
rz(1.705709) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4201323) q[2];
sx q[2];
rz(-0.74154343) q[2];
sx q[2];
rz(0.78765819) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1159462) q[1];
sx q[1];
rz(-1.0315511) q[1];
sx q[1];
rz(2.715452) q[1];
x q[2];
rz(-3.065751) q[3];
sx q[3];
rz(-1.3170529) q[3];
sx q[3];
rz(2.6395519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9289916) q[2];
sx q[2];
rz(-1.5233728) q[2];
sx q[2];
rz(-3.1316481) q[2];
rz(0.11659226) q[3];
sx q[3];
rz(-2.8023585) q[3];
sx q[3];
rz(-0.54178437) q[3];
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
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56124878) q[0];
sx q[0];
rz(-1.9191701) q[0];
sx q[0];
rz(-2.5196581) q[0];
rz(-1.4382582) q[1];
sx q[1];
rz(-0.57241264) q[1];
sx q[1];
rz(0.66351801) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91231649) q[0];
sx q[0];
rz(-1.352542) q[0];
sx q[0];
rz(-2.403232) q[0];
rz(-pi) q[1];
rz(-1.8542669) q[2];
sx q[2];
rz(-1.6890235) q[2];
sx q[2];
rz(0.17519874) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.864526) q[1];
sx q[1];
rz(-0.73523318) q[1];
sx q[1];
rz(-1.0864054) q[1];
rz(-0.31511735) q[3];
sx q[3];
rz(-0.2956008) q[3];
sx q[3];
rz(-2.78418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5180987) q[2];
sx q[2];
rz(-1.3891209) q[2];
sx q[2];
rz(-2.7790879) q[2];
rz(2.6643961) q[3];
sx q[3];
rz(-2.0878744) q[3];
sx q[3];
rz(1.1227192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.057673205) q[0];
sx q[0];
rz(-0.73129439) q[0];
sx q[0];
rz(2.2398563) q[0];
rz(0.67507356) q[1];
sx q[1];
rz(-0.98552862) q[1];
sx q[1];
rz(-2.9973082) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0269723) q[0];
sx q[0];
rz(-1.2204224) q[0];
sx q[0];
rz(-0.087894126) q[0];
rz(-pi) q[1];
rz(-0.4298019) q[2];
sx q[2];
rz(-0.24082213) q[2];
sx q[2];
rz(-0.69320971) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5965138) q[1];
sx q[1];
rz(-0.94280134) q[1];
sx q[1];
rz(1.8819811) q[1];
rz(-1.7606401) q[3];
sx q[3];
rz(-1.1355054) q[3];
sx q[3];
rz(-1.3224885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5054063) q[2];
sx q[2];
rz(-1.212684) q[2];
sx q[2];
rz(-2.9409161) q[2];
rz(-2.0319273) q[3];
sx q[3];
rz(-2.6789013) q[3];
sx q[3];
rz(2.5698575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5526445) q[0];
sx q[0];
rz(-0.968796) q[0];
sx q[0];
rz(2.0198685) q[0];
rz(2.6523392) q[1];
sx q[1];
rz(-1.6901292) q[1];
sx q[1];
rz(2.1314175) q[1];
rz(-0.092838661) q[2];
sx q[2];
rz(-1.9839109) q[2];
sx q[2];
rz(-2.4674923) q[2];
rz(-0.44224593) q[3];
sx q[3];
rz(-1.5968235) q[3];
sx q[3];
rz(2.1376283) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
