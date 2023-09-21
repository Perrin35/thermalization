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
rz(-1.36716258525848) q[0];
sx q[0];
rz(4.05515203078324) q[0];
sx q[0];
rz(11.1542963743131) q[0];
rz(-2.98677587509155) q[1];
sx q[1];
rz(5.68754163582856) q[1];
sx q[1];
rz(7.7653801202695) q[1];
cx q[1],q[0];
rz(0.952811300754547) q[0];
sx q[0];
rz(2.6872098167711) q[0];
sx q[0];
rz(13.7768211126248) q[0];
rz(1.10112321376801) q[2];
sx q[2];
rz(2.85475197632844) q[2];
sx q[2];
rz(4.6341647863309) q[2];
cx q[2],q[1];
rz(-0.239562422037125) q[1];
sx q[1];
rz(8.75034299691255) q[1];
sx q[1];
rz(10.1412954091947) q[1];
rz(0.159266024827957) q[3];
sx q[3];
rz(5.71461644967134) q[3];
sx q[3];
rz(5.14412400721713) q[3];
cx q[3],q[2];
rz(5.29807901382446) q[2];
sx q[2];
rz(0.509223135309764) q[2];
sx q[2];
rz(8.55896071194812) q[2];
rz(-2.18728971481323) q[3];
sx q[3];
rz(4.67998960812623) q[3];
sx q[3];
rz(10.7125643253247) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.998254776000977) q[0];
sx q[0];
rz(4.57823589642579) q[0];
sx q[0];
rz(9.39855862817868) q[0];
rz(1.54013085365295) q[1];
sx q[1];
rz(4.68435052235658) q[1];
sx q[1];
rz(14.7444839239041) q[1];
cx q[1],q[0];
rz(3.16857528686523) q[0];
sx q[0];
rz(6.89722743828828) q[0];
sx q[0];
rz(9.42877056020453) q[0];
rz(1.24047040939331) q[2];
sx q[2];
rz(4.16842988331849) q[2];
sx q[2];
rz(7.04056069850131) q[2];
cx q[2],q[1];
rz(0.023167921230197) q[1];
sx q[1];
rz(4.39284911950166) q[1];
sx q[1];
rz(9.9514956831853) q[1];
rz(-0.601996839046478) q[3];
sx q[3];
rz(0.417998226481028) q[3];
sx q[3];
rz(9.22660153209373) q[3];
cx q[3],q[2];
rz(1.62712371349335) q[2];
sx q[2];
rz(2.01419785817201) q[2];
sx q[2];
rz(12.4318465948026) q[2];
rz(-2.39658045768738) q[3];
sx q[3];
rz(6.5101274569803) q[3];
sx q[3];
rz(7.22591850756809) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.92982518672943) q[0];
sx q[0];
rz(2.75245017011697) q[0];
sx q[0];
rz(7.08062217234775) q[0];
rz(2.09393167495728) q[1];
sx q[1];
rz(6.13344755967195) q[1];
sx q[1];
rz(8.8647812962453) q[1];
cx q[1],q[0];
rz(-3.00217342376709) q[0];
sx q[0];
rz(8.23612943490083) q[0];
sx q[0];
rz(4.49427363871738) q[0];
rz(0.0213039051741362) q[2];
sx q[2];
rz(7.55079761345918) q[2];
sx q[2];
rz(1.49591777323886) q[2];
cx q[2],q[1];
rz(-4.45419645309448) q[1];
sx q[1];
rz(4.51394418080384) q[1];
sx q[1];
rz(8.2203658580701) q[1];
rz(-2.65817260742188) q[3];
sx q[3];
rz(5.1410278399759) q[3];
sx q[3];
rz(9.02407110332652) q[3];
cx q[3],q[2];
rz(2.38931632041931) q[2];
sx q[2];
rz(1.19767776330049) q[2];
sx q[2];
rz(9.5973134547393) q[2];
rz(-0.982073843479156) q[3];
sx q[3];
rz(4.538602860766) q[3];
sx q[3];
rz(10.48274741172) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.13837063312531) q[0];
sx q[0];
rz(5.18556681473786) q[0];
sx q[0];
rz(9.70929229854747) q[0];
rz(6.59988927841187) q[1];
sx q[1];
rz(5.85042205651338) q[1];
sx q[1];
rz(4.44032380580112) q[1];
cx q[1],q[0];
rz(0.117677800357342) q[0];
sx q[0];
rz(2.06645837624604) q[0];
sx q[0];
rz(10.7390690803449) q[0];
rz(-3.90800380706787) q[2];
sx q[2];
rz(1.57172361214692) q[2];
sx q[2];
rz(4.73317334651157) q[2];
cx q[2],q[1];
rz(2.27007699012756) q[1];
sx q[1];
rz(4.5260752757364) q[1];
sx q[1];
rz(8.85286668538257) q[1];
rz(0.32818016409874) q[3];
sx q[3];
rz(1.67718520958955) q[3];
sx q[3];
rz(8.42193696498081) q[3];
cx q[3],q[2];
rz(1.60568845272064) q[2];
sx q[2];
rz(3.33742609818513) q[2];
sx q[2];
rz(5.8985051870267) q[2];
rz(3.89565205574036) q[3];
sx q[3];
rz(4.19810965855653) q[3];
sx q[3];
rz(7.73748741149112) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.53831791877747) q[0];
sx q[0];
rz(4.12869259913499) q[0];
sx q[0];
rz(11.1797144174497) q[0];
rz(0.231001347303391) q[1];
sx q[1];
rz(4.94203129609162) q[1];
sx q[1];
rz(9.72160456179782) q[1];
cx q[1],q[0];
rz(3.55403256416321) q[0];
sx q[0];
rz(4.67322567303712) q[0];
sx q[0];
rz(8.42504600285693) q[0];
rz(0.827956736087799) q[2];
sx q[2];
rz(4.82546368439729) q[2];
sx q[2];
rz(9.8856485247533) q[2];
cx q[2],q[1];
rz(1.24248027801514) q[1];
sx q[1];
rz(5.65092268784577) q[1];
sx q[1];
rz(8.71507195233508) q[1];
rz(0.0889256224036217) q[3];
sx q[3];
rz(4.00067177613313) q[3];
sx q[3];
rz(12.0199427366178) q[3];
cx q[3],q[2];
rz(2.76328730583191) q[2];
sx q[2];
rz(1.83132353623445) q[2];
sx q[2];
rz(5.89070937632724) q[2];
rz(-1.15224194526672) q[3];
sx q[3];
rz(3.85617986519868) q[3];
sx q[3];
rz(12.8838150262754) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.015793800354) q[0];
sx q[0];
rz(1.56904652913148) q[0];
sx q[0];
rz(8.67339197396442) q[0];
rz(-1.32793509960175) q[1];
sx q[1];
rz(4.40498057206208) q[1];
sx q[1];
rz(10.031117117397) q[1];
cx q[1],q[0];
rz(-0.408892601728439) q[0];
sx q[0];
rz(3.09282867808873) q[0];
sx q[0];
rz(11.3439632415692) q[0];
rz(-1.15433204174042) q[2];
sx q[2];
rz(0.935668619471141) q[2];
sx q[2];
rz(6.37658808230563) q[2];
cx q[2],q[1];
rz(2.63410830497742) q[1];
sx q[1];
rz(1.02352300484712) q[1];
sx q[1];
rz(7.64230010508701) q[1];
rz(1.64579617977142) q[3];
sx q[3];
rz(4.92936769326264) q[3];
sx q[3];
rz(7.33554408549472) q[3];
cx q[3],q[2];
rz(-0.638858258724213) q[2];
sx q[2];
rz(7.32493177254731) q[2];
sx q[2];
rz(7.4225222825925) q[2];
rz(-1.48494386672974) q[3];
sx q[3];
rz(1.96102014382417) q[3];
sx q[3];
rz(9.32052073477908) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.60956025123596) q[0];
sx q[0];
rz(0.748152645426341) q[0];
sx q[0];
rz(8.91667739152118) q[0];
rz(-1.56289005279541) q[1];
sx q[1];
rz(7.37201658089692) q[1];
sx q[1];
rz(10.2150261759679) q[1];
cx q[1],q[0];
rz(3.3398711681366) q[0];
sx q[0];
rz(0.675307901697703) q[0];
sx q[0];
rz(8.33792803286716) q[0];
rz(1.51677072048187) q[2];
sx q[2];
rz(6.49712887604768) q[2];
sx q[2];
rz(21.6592607259671) q[2];
cx q[2],q[1];
rz(2.02231693267822) q[1];
sx q[1];
rz(4.32225445111329) q[1];
sx q[1];
rz(12.2739608049314) q[1];
rz(3.51924133300781) q[3];
sx q[3];
rz(2.18340632517869) q[3];
sx q[3];
rz(8.53713575600787) q[3];
cx q[3],q[2];
rz(3.19467782974243) q[2];
sx q[2];
rz(3.58343967993791) q[2];
sx q[2];
rz(10.8379947900693) q[2];
rz(-1.66480708122253) q[3];
sx q[3];
rz(1.03527417977388) q[3];
sx q[3];
rz(12.6279251336972) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.42478942871094) q[0];
sx q[0];
rz(3.17190348555381) q[0];
sx q[0];
rz(10.4720042705457) q[0];
rz(0.609102427959442) q[1];
sx q[1];
rz(4.86921647389466) q[1];
sx q[1];
rz(11.1776579379956) q[1];
cx q[1],q[0];
rz(3.26185655593872) q[0];
sx q[0];
rz(4.74219587643678) q[0];
sx q[0];
rz(12.006735777847) q[0];
rz(7.19532537460327) q[2];
sx q[2];
rz(5.64687887032563) q[2];
sx q[2];
rz(12.6181106328885) q[2];
cx q[2],q[1];
rz(0.705062389373779) q[1];
sx q[1];
rz(0.0925180037790021) q[1];
sx q[1];
rz(7.73047683238193) q[1];
rz(2.72693610191345) q[3];
sx q[3];
rz(3.78998097975785) q[3];
sx q[3];
rz(8.78875843285724) q[3];
cx q[3],q[2];
rz(-1.95288145542145) q[2];
sx q[2];
rz(5.87295022805268) q[2];
sx q[2];
rz(8.54252258538409) q[2];
rz(-1.40117180347443) q[3];
sx q[3];
rz(4.30795029004151) q[3];
sx q[3];
rz(8.22428319453403) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.271544873714447) q[0];
sx q[0];
rz(5.86637941201264) q[0];
sx q[0];
rz(7.99868414401218) q[0];
rz(3.22305393218994) q[1];
sx q[1];
rz(4.30416086514527) q[1];
sx q[1];
rz(8.86654517649814) q[1];
cx q[1],q[0];
rz(-0.249606251716614) q[0];
sx q[0];
rz(5.68925729592378) q[0];
sx q[0];
rz(10.165398156635) q[0];
rz(-3.24906086921692) q[2];
sx q[2];
rz(1.62459138234193) q[2];
sx q[2];
rz(7.97981367110416) q[2];
cx q[2],q[1];
rz(-3.26962208747864) q[1];
sx q[1];
rz(5.45470705826814) q[1];
sx q[1];
rz(14.2238712072293) q[1];
rz(-1.22983598709106) q[3];
sx q[3];
rz(3.65294936497743) q[3];
sx q[3];
rz(11.6627838373105) q[3];
cx q[3],q[2];
rz(1.84908628463745) q[2];
sx q[2];
rz(4.41092005570466) q[2];
sx q[2];
rz(9.21259396373435) q[2];
rz(2.92961812019348) q[3];
sx q[3];
rz(3.82484981616075) q[3];
sx q[3];
rz(10.6268262624662) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.504871726036072) q[0];
sx q[0];
rz(2.31754055817659) q[0];
sx q[0];
rz(11.028531408302) q[0];
rz(-3.96699976921082) q[1];
sx q[1];
rz(5.61041918595368) q[1];
sx q[1];
rz(6.80647990702792) q[1];
cx q[1],q[0];
rz(2.97420454025269) q[0];
sx q[0];
rz(3.75358841021592) q[0];
sx q[0];
rz(7.96228883265659) q[0];
rz(-1.5269581079483) q[2];
sx q[2];
rz(5.28390422661836) q[2];
sx q[2];
rz(9.72063515185519) q[2];
cx q[2],q[1];
rz(1.11124217510223) q[1];
sx q[1];
rz(5.79673901398713) q[1];
sx q[1];
rz(10.2727763414304) q[1];
rz(-1.81577682495117) q[3];
sx q[3];
rz(1.79532006581361) q[3];
sx q[3];
rz(8.97372916936084) q[3];
cx q[3],q[2];
rz(-1.30451965332031) q[2];
sx q[2];
rz(3.85707375605638) q[2];
sx q[2];
rz(9.69408095478221) q[2];
rz(-0.494201600551605) q[3];
sx q[3];
rz(3.98794958193833) q[3];
sx q[3];
rz(11.4802808523099) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.221064358949661) q[0];
sx q[0];
rz(6.19277349312837) q[0];
sx q[0];
rz(8.95699990390941) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(1.12656760215759) q[1];
sx q[1];
rz(1.2943892796808) q[1];
sx q[1];
rz(7.75759909152194) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(0.492837995290756) q[2];
sx q[2];
rz(1.34314122994477) q[2];
sx q[2];
rz(7.98679778575107) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(1.07378196716309) q[3];
sx q[3];
rz(0.976461323099681) q[3];
sx q[3];
rz(8.40702483653232) q[3];
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