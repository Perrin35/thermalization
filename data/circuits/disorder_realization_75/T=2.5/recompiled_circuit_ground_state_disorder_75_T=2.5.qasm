OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg meas[4];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(pi/2) q[3];
sx q[3];
rz(3*pi/2) q[3];
sx q[3];
rz(5*pi/2) q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.925754427909851) q[0];
sx q[0];
rz(2.8211419900232) q[0];
sx q[0];
rz(9.55316432415649) q[0];
rz(4.09659624099731) q[1];
sx q[1];
rz(8.76488653023774) q[1];
sx q[1];
rz(13.1534912347715) q[1];
cx q[1],q[0];
rz(-0.05953723564744) q[0];
sx q[0];
rz(4.24794426758821) q[0];
sx q[0];
rz(12.4046995401303) q[0];
rz(-5.41933012008667) q[2];
sx q[2];
rz(3.79078450997407) q[2];
sx q[2];
rz(9.87714884280368) q[2];
cx q[2],q[1];
rz(0.942266643047333) q[1];
sx q[1];
rz(5.14713421662385) q[1];
sx q[1];
rz(8.49600855111285) q[1];
rz(-3.13495659828186) q[3];
sx q[3];
rz(9.0402833541208) q[3];
sx q[3];
rz(15.8746576070706) q[3];
cx q[3],q[2];
rz(6.54855060577393) q[2];
sx q[2];
rz(2.90761248965795) q[2];
sx q[2];
rz(8.50375286339923) q[2];
rz(4.66625022888184) q[3];
sx q[3];
rz(4.99205950100953) q[3];
sx q[3];
rz(6.0974638223569) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.281445354223251) q[0];
sx q[0];
rz(3.53650337656076) q[0];
sx q[0];
rz(10.7741559505384) q[0];
rz(-0.347406268119812) q[1];
sx q[1];
rz(5.09432307084138) q[1];
sx q[1];
rz(10.8278171777646) q[1];
cx q[1],q[0];
rz(-3.27737784385681) q[0];
sx q[0];
rz(6.91135445435578) q[0];
sx q[0];
rz(13.7632136106412) q[0];
rz(-0.61739057302475) q[2];
sx q[2];
rz(7.13298765023286) q[2];
sx q[2];
rz(4.55796525477573) q[2];
cx q[2],q[1];
rz(2.78293204307556) q[1];
sx q[1];
rz(0.730507763224193) q[1];
sx q[1];
rz(10.7647356748502) q[1];
rz(-2.81770968437195) q[3];
sx q[3];
rz(6.99773231347138) q[3];
sx q[3];
rz(5.16319749354526) q[3];
cx q[3],q[2];
rz(-1.81895244121552) q[2];
sx q[2];
rz(4.44945839245851) q[2];
sx q[2];
rz(9.5835999905984) q[2];
rz(-3.35924768447876) q[3];
sx q[3];
rz(1.75261512597138) q[3];
sx q[3];
rz(7.63369605540439) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.366671025753021) q[0];
sx q[0];
rz(-1.3008302132315) q[0];
sx q[0];
rz(11.2969033479612) q[0];
rz(-0.419957548379898) q[1];
sx q[1];
rz(5.28234401543672) q[1];
sx q[1];
rz(4.74397704600498) q[1];
cx q[1],q[0];
rz(3.54678821563721) q[0];
sx q[0];
rz(4.27087715466554) q[0];
sx q[0];
rz(13.7910895109098) q[0];
rz(-3.32526206970215) q[2];
sx q[2];
rz(1.35461262066896) q[2];
sx q[2];
rz(8.77619496583148) q[2];
cx q[2],q[1];
rz(-0.532454133033752) q[1];
sx q[1];
rz(-1.62416061560576) q[1];
sx q[1];
rz(6.67878911494418) q[1];
rz(-1.07755875587463) q[3];
sx q[3];
rz(6.90892687638337) q[3];
sx q[3];
rz(10.8566030025403) q[3];
cx q[3],q[2];
rz(0.0821082144975662) q[2];
sx q[2];
rz(2.29936453898484) q[2];
sx q[2];
rz(12.3746716737668) q[2];
rz(2.58904027938843) q[3];
sx q[3];
rz(4.66660919983918) q[3];
sx q[3];
rz(11.54189488887) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.73156523704529) q[0];
sx q[0];
rz(6.93160858948762) q[0];
sx q[0];
rz(6.88755056857272) q[0];
rz(1.89616894721985) q[1];
sx q[1];
rz(5.53289237816865) q[1];
sx q[1];
rz(7.1429321527402) q[1];
cx q[1],q[0];
rz(-2.51675462722778) q[0];
sx q[0];
rz(5.55838385422761) q[0];
sx q[0];
rz(9.94764325617954) q[0];
rz(3.0395200252533) q[2];
sx q[2];
rz(1.58491757710511) q[2];
sx q[2];
rz(9.95700196026965) q[2];
cx q[2],q[1];
rz(3.68956208229065) q[1];
sx q[1];
rz(7.67242017586763) q[1];
sx q[1];
rz(9.47538182734653) q[1];
rz(2.88799738883972) q[3];
sx q[3];
rz(5.41665664513642) q[3];
sx q[3];
rz(9.63715634345218) q[3];
cx q[3],q[2];
rz(1.8613646030426) q[2];
sx q[2];
rz(2.12755289872224) q[2];
sx q[2];
rz(7.72178361415073) q[2];
rz(1.29226183891296) q[3];
sx q[3];
rz(3.96946403582627) q[3];
sx q[3];
rz(4.19621512889072) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.07865214347839) q[0];
sx q[0];
rz(3.41473990877206) q[0];
sx q[0];
rz(9.10692737101718) q[0];
rz(1.33927965164185) q[1];
sx q[1];
rz(-2.72367176215117) q[1];
sx q[1];
rz(8.44606182574435) q[1];
cx q[1],q[0];
rz(-0.624171078205109) q[0];
sx q[0];
rz(4.68213561375672) q[0];
sx q[0];
rz(7.45291683673068) q[0];
rz(-1.37543618679047) q[2];
sx q[2];
rz(1.4835411628061) q[2];
sx q[2];
rz(10.7877919435422) q[2];
cx q[2],q[1];
rz(1.38263440132141) q[1];
sx q[1];
rz(0.823732288675853) q[1];
sx q[1];
rz(9.252751967303) q[1];
rz(-2.35424876213074) q[3];
sx q[3];
rz(3.41132861574227) q[3];
sx q[3];
rz(8.07774958609744) q[3];
cx q[3],q[2];
rz(6.481529712677) q[2];
sx q[2];
rz(5.33815422852571) q[2];
sx q[2];
rz(16.8252601385038) q[2];
rz(-0.18103389441967) q[3];
sx q[3];
rz(5.04625371296937) q[3];
sx q[3];
rz(7.48044905661746) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-3.84659790992737) q[0];
sx q[0];
rz(6.01610341866548) q[0];
sx q[0];
rz(11.2566891670148) q[0];
rz(-7.30505561828613) q[1];
sx q[1];
rz(5.44649830658967) q[1];
sx q[1];
rz(7.77172157763644) q[1];
cx q[1],q[0];
rz(-0.903279006481171) q[0];
sx q[0];
rz(5.36146298249299) q[0];
sx q[0];
rz(11.4179693221967) q[0];
rz(-4.99953746795654) q[2];
sx q[2];
rz(3.84040811856324) q[2];
sx q[2];
rz(12.4407899141233) q[2];
cx q[2],q[1];
rz(5.75636577606201) q[1];
sx q[1];
rz(2.0627910216623) q[1];
sx q[1];
rz(9.58658654092952) q[1];
rz(-1.68453121185303) q[3];
sx q[3];
rz(4.1544430573755) q[3];
sx q[3];
rz(6.5008389711301) q[3];
cx q[3],q[2];
rz(-3.80429458618164) q[2];
sx q[2];
rz(1.75444892247254) q[2];
sx q[2];
rz(12.3022889852445) q[2];
rz(-4.61798334121704) q[3];
sx q[3];
rz(5.6030332167917) q[3];
sx q[3];
rz(9.8445852458398) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.611780941486359) q[0];
sx q[0];
rz(4.78229955037171) q[0];
sx q[0];
rz(10.4890379667203) q[0];
rz(3.24078941345215) q[1];
sx q[1];
rz(1.79250374634797) q[1];
sx q[1];
rz(7.96420631407901) q[1];
cx q[1],q[0];
rz(1.19535803794861) q[0];
sx q[0];
rz(4.77615884144837) q[0];
sx q[0];
rz(7.36920545100375) q[0];
rz(0.12701290845871) q[2];
sx q[2];
rz(3.44268232782418) q[2];
sx q[2];
rz(4.95951554774448) q[2];
cx q[2],q[1];
rz(3.48173451423645) q[1];
sx q[1];
rz(3.65507349570329) q[1];
sx q[1];
rz(7.25047371386691) q[1];
rz(2.25494170188904) q[3];
sx q[3];
rz(5.64105621178681) q[3];
sx q[3];
rz(7.0475337266843) q[3];
cx q[3],q[2];
rz(1.35082411766052) q[2];
sx q[2];
rz(2.02275732358033) q[2];
sx q[2];
rz(11.9668724298398) q[2];
rz(-5.38553953170776) q[3];
sx q[3];
rz(4.86362531979615) q[3];
sx q[3];
rz(9.38676310553356) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.813832342624664) q[0];
sx q[0];
rz(4.20194378693635) q[0];
sx q[0];
rz(8.02333376406833) q[0];
rz(0.0714881494641304) q[1];
sx q[1];
rz(4.60495260556275) q[1];
sx q[1];
rz(7.28432533740207) q[1];
cx q[1],q[0];
rz(-1.63426959514618) q[0];
sx q[0];
rz(1.8232913335138) q[0];
sx q[0];
rz(13.6758055448453) q[0];
rz(2.75433158874512) q[2];
sx q[2];
rz(5.33908048470552) q[2];
sx q[2];
rz(10.7284996271054) q[2];
cx q[2],q[1];
rz(2.04671716690063) q[1];
sx q[1];
rz(6.93015805085237) q[1];
sx q[1];
rz(13.9173717260282) q[1];
rz(0.866451561450958) q[3];
sx q[3];
rz(2.33726903994615) q[3];
sx q[3];
rz(11.8104083299558) q[3];
cx q[3],q[2];
rz(-5.1059684753418) q[2];
sx q[2];
rz(1.55286768277223) q[2];
sx q[2];
rz(13.2848114728849) q[2];
rz(0.0467723794281483) q[3];
sx q[3];
rz(-0.707435695332936) q[3];
sx q[3];
rz(8.85497990845844) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.30963322520256) q[0];
sx q[0];
rz(3.29600408871705) q[0];
sx q[0];
rz(13.256818509094) q[0];
rz(0.180282324552536) q[1];
sx q[1];
rz(4.20281914074952) q[1];
sx q[1];
rz(5.96043512820407) q[1];
cx q[1],q[0];
rz(-2.66473412513733) q[0];
sx q[0];
rz(3.32799191971356) q[0];
sx q[0];
rz(7.25295708178684) q[0];
rz(2.58947539329529) q[2];
sx q[2];
rz(4.7068796475702) q[2];
sx q[2];
rz(7.95111320017978) q[2];
cx q[2],q[1];
rz(-1.91959357261658) q[1];
sx q[1];
rz(1.02424350579316) q[1];
sx q[1];
rz(11.9291212320249) q[1];
rz(0.812028169631958) q[3];
sx q[3];
rz(3.98226365645463) q[3];
sx q[3];
rz(7.44278297423526) q[3];
cx q[3],q[2];
rz(3.61260843276978) q[2];
sx q[2];
rz(-1.34795030753081) q[2];
sx q[2];
rz(13.5842523336331) q[2];
rz(-0.00353035912849009) q[3];
sx q[3];
rz(4.11120197375352) q[3];
sx q[3];
rz(11.3544790506284) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.328053504228592) q[0];
sx q[0];
rz(0.737547548609324) q[0];
sx q[0];
rz(10.3507122754972) q[0];
rz(-0.137358590960503) q[1];
sx q[1];
rz(0.672471435862132) q[1];
sx q[1];
rz(5.68327352999851) q[1];
cx q[1],q[0];
rz(-3.16533541679382) q[0];
sx q[0];
rz(4.70968464215333) q[0];
sx q[0];
rz(10.4516440391461) q[0];
rz(2.27657222747803) q[2];
sx q[2];
rz(3.90453621943528) q[2];
sx q[2];
rz(8.3145642042081) q[2];
cx q[2],q[1];
rz(-2.74739003181458) q[1];
sx q[1];
rz(8.83530393441255) q[1];
sx q[1];
rz(12.5116193056028) q[1];
rz(0.581688642501831) q[3];
sx q[3];
rz(5.73264876206452) q[3];
sx q[3];
rz(10.8607474327008) q[3];
cx q[3],q[2];
rz(7.69565868377686) q[2];
sx q[2];
rz(5.26425877411897) q[2];
sx q[2];
rz(2.49983832835361) q[2];
rz(1.98520910739899) q[3];
sx q[3];
rz(3.90280059178407) q[3];
sx q[3];
rz(7.80644927024051) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-2.20875382423401) q[0];
sx q[0];
rz(2.14090600808198) q[0];
sx q[0];
rz(10.2836165189664) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(0.137821272015572) q[1];
sx q[1];
rz(7.78548303444917) q[1];
sx q[1];
rz(6.38275454043552) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(1.52957129478455) q[2];
sx q[2];
rz(1.81943717797334) q[2];
sx q[2];
rz(11.1986917018811) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(-3.29929876327515) q[3];
sx q[3];
rz(2.17443028290803) q[3];
sx q[3];
rz(7.50673840045139) q[3];
rz(pi/2) q[3];
sx q[3];
rz(3*pi/2) q[3];
sx q[3];
rz(5*pi/2) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> meas[0];
measure q[1] -> meas[1];
measure q[2] -> meas[2];
measure q[3] -> meas[3];
