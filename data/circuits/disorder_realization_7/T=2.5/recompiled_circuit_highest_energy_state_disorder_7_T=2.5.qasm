OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg meas[4];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
sx q[3];
rz(5*pi/2) q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.46472907066345) q[0];
sx q[0];
rz(1.61141768296296) q[0];
sx q[0];
rz(9.72206977605029) q[0];
rz(-5.50414323806763) q[1];
sx q[1];
rz(6.12255230744416) q[1];
sx q[1];
rz(10.8607220411222) q[1];
cx q[1],q[0];
rz(2.44776177406311) q[0];
sx q[0];
rz(4.20743826230104) q[0];
sx q[0];
rz(8.79116765259906) q[0];
rz(-1.4421671628952) q[2];
sx q[2];
rz(1.82782152493531) q[2];
sx q[2];
rz(10.2687653064649) q[2];
cx q[2],q[1];
rz(-1.04321920871735) q[1];
sx q[1];
rz(2.34651521046693) q[1];
sx q[1];
rz(14.5951204061429) q[1];
rz(4.14607048034668) q[3];
sx q[3];
rz(5.92227497895295) q[3];
sx q[3];
rz(13.2147817373197) q[3];
cx q[3],q[2];
rz(1.54331886768341) q[2];
sx q[2];
rz(2.71365493734414) q[2];
sx q[2];
rz(9.59256804584667) q[2];
rz(1.12818908691406) q[3];
sx q[3];
rz(4.71066096623475) q[3];
sx q[3];
rz(11.3771620750348) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.0968649163842201) q[0];
sx q[0];
rz(3.81172755559022) q[0];
sx q[0];
rz(8.34973857402011) q[0];
rz(-1.77544045448303) q[1];
sx q[1];
rz(1.75518825848634) q[1];
sx q[1];
rz(10.7385600566785) q[1];
cx q[1],q[0];
rz(0.157206505537033) q[0];
sx q[0];
rz(4.4909304698282) q[0];
sx q[0];
rz(10.6588488578717) q[0];
rz(1.38448715209961) q[2];
sx q[2];
rz(4.54583779175813) q[2];
sx q[2];
rz(9.32973705082341) q[2];
cx q[2],q[1];
rz(-2.94568538665771) q[1];
sx q[1];
rz(5.43422356446321) q[1];
sx q[1];
rz(14.0348491430204) q[1];
rz(1.70442080497742) q[3];
sx q[3];
rz(5.04058984120423) q[3];
sx q[3];
rz(11.092083311073) q[3];
cx q[3],q[2];
rz(3.79269862174988) q[2];
sx q[2];
rz(7.65969672997529) q[2];
sx q[2];
rz(9.00039309858485) q[2];
rz(-0.732447981834412) q[3];
sx q[3];
rz(4.6248278935724) q[3];
sx q[3];
rz(10.8881939411084) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.803777933120728) q[0];
sx q[0];
rz(3.85659536917741) q[0];
sx q[0];
rz(13.0691404104154) q[0];
rz(-4.08523464202881) q[1];
sx q[1];
rz(3.7840452512079) q[1];
sx q[1];
rz(11.7511772870938) q[1];
cx q[1],q[0];
rz(1.07499325275421) q[0];
sx q[0];
rz(2.94783589442308) q[0];
sx q[0];
rz(10.7501293182294) q[0];
rz(3.9853675365448) q[2];
sx q[2];
rz(2.98248483439023) q[2];
sx q[2];
rz(6.95010206698581) q[2];
cx q[2],q[1];
rz(0.803226411342621) q[1];
sx q[1];
rz(5.91002074082429) q[1];
sx q[1];
rz(10.205433702461) q[1];
rz(1.00913095474243) q[3];
sx q[3];
rz(5.32917800744111) q[3];
sx q[3];
rz(7.63464710711643) q[3];
cx q[3],q[2];
rz(-0.0766855850815773) q[2];
sx q[2];
rz(4.83205238183076) q[2];
sx q[2];
rz(8.29827044009372) q[2];
rz(-1.35697758197784) q[3];
sx q[3];
rz(6.78280440171296) q[3];
sx q[3];
rz(12.4264321088712) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.559468805789948) q[0];
sx q[0];
rz(2.42156502802903) q[0];
sx q[0];
rz(9.06486404537364) q[0];
rz(3.38238072395325) q[1];
sx q[1];
rz(1.99379316170747) q[1];
sx q[1];
rz(9.05158857106372) q[1];
cx q[1],q[0];
rz(0.686158955097198) q[0];
sx q[0];
rz(2.96692744095857) q[0];
sx q[0];
rz(10.0311794042508) q[0];
rz(0.234747469425201) q[2];
sx q[2];
rz(2.12146070797975) q[2];
sx q[2];
rz(13.1178819894712) q[2];
cx q[2],q[1];
rz(-3.28387522697449) q[1];
sx q[1];
rz(5.0526010115915) q[1];
sx q[1];
rz(11.545834517471) q[1];
rz(2.95110654830933) q[3];
sx q[3];
rz(-1.13087400595611) q[3];
sx q[3];
rz(11.6751856565396) q[3];
cx q[3],q[2];
rz(-0.0844262465834618) q[2];
sx q[2];
rz(4.57361486752565) q[2];
sx q[2];
rz(8.51633474826022) q[2];
rz(4.21157169342041) q[3];
sx q[3];
rz(5.09862020810182) q[3];
sx q[3];
rz(7.26647398471042) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.958036065101624) q[0];
sx q[0];
rz(4.03520110447938) q[0];
sx q[0];
rz(8.51173368691608) q[0];
rz(-2.45195817947388) q[1];
sx q[1];
rz(5.37441125710542) q[1];
sx q[1];
rz(8.65076617001697) q[1];
cx q[1],q[0];
rz(0.900923132896423) q[0];
sx q[0];
rz(4.25604060490663) q[0];
sx q[0];
rz(7.19293353556796) q[0];
rz(-3.10985660552979) q[2];
sx q[2];
rz(5.50772419770295) q[2];
sx q[2];
rz(9.68547282218143) q[2];
cx q[2],q[1];
rz(1.41786313056946) q[1];
sx q[1];
rz(3.61501339276368) q[1];
sx q[1];
rz(9.85406929849788) q[1];
rz(-1.02294862270355) q[3];
sx q[3];
rz(1.59611478646333) q[3];
sx q[3];
rz(11.930621123306) q[3];
cx q[3],q[2];
rz(-1.24560499191284) q[2];
sx q[2];
rz(1.39348033268983) q[2];
sx q[2];
rz(8.68527863024875) q[2];
rz(1.63515222072601) q[3];
sx q[3];
rz(2.20619008143479) q[3];
sx q[3];
rz(8.61705020665332) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.62261152267456) q[0];
sx q[0];
rz(3.93672815163667) q[0];
sx q[0];
rz(8.29948899745151) q[0];
rz(3.00334239006042) q[1];
sx q[1];
rz(4.81592026551301) q[1];
sx q[1];
rz(7.85757241248294) q[1];
cx q[1],q[0];
rz(1.34309053421021) q[0];
sx q[0];
rz(6.78416934807832) q[0];
sx q[0];
rz(13.4083759546201) q[0];
rz(2.84944796562195) q[2];
sx q[2];
rz(5.110140593844) q[2];
sx q[2];
rz(11.0393150806348) q[2];
cx q[2],q[1];
rz(4.38281965255737) q[1];
sx q[1];
rz(4.45396688778932) q[1];
sx q[1];
rz(12.4377548456113) q[1];
rz(-1.02155065536499) q[3];
sx q[3];
rz(4.54327085812623) q[3];
sx q[3];
rz(8.36453721522495) q[3];
cx q[3],q[2];
rz(-0.773419976234436) q[2];
sx q[2];
rz(8.3566953261667) q[2];
sx q[2];
rz(14.6087727308194) q[2];
rz(-2.15589714050293) q[3];
sx q[3];
rz(4.76699033577973) q[3];
sx q[3];
rz(12.1337859392087) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.84283351898193) q[0];
sx q[0];
rz(4.13660404284532) q[0];
sx q[0];
rz(13.3751697301786) q[0];
rz(2.47592997550964) q[1];
sx q[1];
rz(5.30368915398652) q[1];
sx q[1];
rz(8.44336549042865) q[1];
cx q[1],q[0];
rz(-0.279413759708405) q[0];
sx q[0];
rz(5.28272500832612) q[0];
sx q[0];
rz(7.92191359996005) q[0];
rz(1.12923622131348) q[2];
sx q[2];
rz(0.549166353540965) q[2];
sx q[2];
rz(14.1205687284391) q[2];
cx q[2],q[1];
rz(-1.04048240184784) q[1];
sx q[1];
rz(5.79698601563508) q[1];
sx q[1];
rz(6.57683775424167) q[1];
rz(3.31336641311646) q[3];
sx q[3];
rz(2.46187588770921) q[3];
sx q[3];
rz(7.62417123316928) q[3];
cx q[3],q[2];
rz(-1.07706463336945) q[2];
sx q[2];
rz(1.47461024125154) q[2];
sx q[2];
rz(9.05576056837245) q[2];
rz(-1.40242326259613) q[3];
sx q[3];
rz(4.02874323924119) q[3];
sx q[3];
rz(10.7460254192273) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.30826652050018) q[0];
sx q[0];
rz(5.45431176026399) q[0];
sx q[0];
rz(9.7549582183282) q[0];
rz(-0.456868320703506) q[1];
sx q[1];
rz(3.57691708405549) q[1];
sx q[1];
rz(10.0230636954229) q[1];
cx q[1],q[0];
rz(-2.90151810646057) q[0];
sx q[0];
rz(5.4759580214792) q[0];
sx q[0];
rz(11.6190573930661) q[0];
rz(-0.469173461198807) q[2];
sx q[2];
rz(4.77568033536012) q[2];
sx q[2];
rz(11.8958191633145) q[2];
cx q[2],q[1];
rz(0.484093248844147) q[1];
sx q[1];
rz(1.34170213540132) q[1];
sx q[1];
rz(11.3377709150235) q[1];
rz(2.55374431610107) q[3];
sx q[3];
rz(1.00252571900422) q[3];
sx q[3];
rz(10.4070818185727) q[3];
cx q[3],q[2];
rz(-0.869970798492432) q[2];
sx q[2];
rz(7.15965238411958) q[2];
sx q[2];
rz(9.48939989357396) q[2];
rz(1.4501074552536) q[3];
sx q[3];
rz(5.88299957116181) q[3];
sx q[3];
rz(14.1120347738187) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.862186312675476) q[0];
sx q[0];
rz(3.47570091684396) q[0];
sx q[0];
rz(8.90351328849002) q[0];
rz(2.08172559738159) q[1];
sx q[1];
rz(1.16185489495332) q[1];
sx q[1];
rz(10.4561567068021) q[1];
cx q[1],q[0];
rz(0.176126211881638) q[0];
sx q[0];
rz(2.30760076840455) q[0];
sx q[0];
rz(12.6665949582975) q[0];
rz(-1.27792227268219) q[2];
sx q[2];
rz(1.68568626244599) q[2];
sx q[2];
rz(9.44889386220976) q[2];
cx q[2],q[1];
rz(-3.18678617477417) q[1];
sx q[1];
rz(0.611224325495311) q[1];
sx q[1];
rz(16.226188158981) q[1];
rz(1.38217902183533) q[3];
sx q[3];
rz(0.484357031183787) q[3];
sx q[3];
rz(12.8175825834195) q[3];
cx q[3],q[2];
rz(1.39702677726746) q[2];
sx q[2];
rz(2.41377452214295) q[2];
sx q[2];
rz(9.33740593343183) q[2];
rz(4.35737895965576) q[3];
sx q[3];
rz(4.60173502762849) q[3];
sx q[3];
rz(12.3962614297788) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.666072964668274) q[0];
sx q[0];
rz(2.25057944853837) q[0];
sx q[0];
rz(8.28561589717075) q[0];
rz(-0.49925222992897) q[1];
sx q[1];
rz(4.59078589280183) q[1];
sx q[1];
rz(11.2329015493314) q[1];
cx q[1],q[0];
rz(0.0697869285941124) q[0];
sx q[0];
rz(3.00203679700429) q[0];
sx q[0];
rz(8.5584572315137) q[0];
rz(-1.52099454402924) q[2];
sx q[2];
rz(4.96703818638856) q[2];
sx q[2];
rz(10.8580204009931) q[2];
cx q[2],q[1];
rz(1.00296103954315) q[1];
sx q[1];
rz(1.31827822526033) q[1];
sx q[1];
rz(9.20213404893085) q[1];
rz(2.31283926963806) q[3];
sx q[3];
rz(2.29314056237275) q[3];
sx q[3];
rz(5.74919674395725) q[3];
cx q[3],q[2];
rz(4.64259767532349) q[2];
sx q[2];
rz(1.86684003670747) q[2];
sx q[2];
rz(6.56642363070651) q[2];
rz(-1.9404274225235) q[3];
sx q[3];
rz(2.46565064986283) q[3];
sx q[3];
rz(10.4702098131101) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.214780002832413) q[0];
sx q[0];
rz(4.03546795447404) q[0];
sx q[0];
rz(11.4112384080808) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(-1.71773433685303) q[1];
sx q[1];
rz(5.8583494742685) q[1];
sx q[1];
rz(9.37810162677571) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(-2.46218633651733) q[2];
sx q[2];
rz(3.69117012818391) q[2];
sx q[2];
rz(10.5761574268262) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(3.11951899528503) q[3];
sx q[3];
rz(3.63917476137216) q[3];
sx q[3];
rz(7.03748438357517) q[3];
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
