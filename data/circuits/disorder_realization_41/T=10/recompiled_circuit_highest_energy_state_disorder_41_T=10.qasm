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
rz(0.323494106531143) q[0];
sx q[0];
rz(2.9383724202686) q[0];
sx q[0];
rz(9.12462244033023) q[0];
rz(-0.35429859161377) q[1];
sx q[1];
rz(4.37591818173463) q[1];
sx q[1];
rz(10.1959225892942) q[1];
cx q[1],q[0];
rz(0.0104926517233253) q[0];
sx q[0];
rz(3.40988186200196) q[0];
sx q[0];
rz(10.218906736366) q[0];
rz(1.03729343414307) q[2];
sx q[2];
rz(4.29637494881684) q[2];
sx q[2];
rz(8.85691366194888) q[2];
cx q[2],q[1];
rz(0.315095454454422) q[1];
sx q[1];
rz(2.57194480498368) q[1];
sx q[1];
rz(10.1704336166303) q[1];
rz(-0.145017340779305) q[3];
sx q[3];
rz(3.57961890299852) q[3];
sx q[3];
rz(11.0407321214597) q[3];
cx q[3],q[2];
rz(0.990420877933502) q[2];
sx q[2];
rz(3.60310349066789) q[2];
sx q[2];
rz(10.5085245132367) q[2];
rz(0.574619650840759) q[3];
sx q[3];
rz(4.73663309414918) q[3];
sx q[3];
rz(10.2021981835286) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.827117800712585) q[0];
sx q[0];
rz(3.7706460078531) q[0];
sx q[0];
rz(9.02629125713512) q[0];
rz(1.99729478359222) q[1];
sx q[1];
rz(2.53563478787477) q[1];
sx q[1];
rz(8.47044900654956) q[1];
cx q[1],q[0];
rz(1.07145571708679) q[0];
sx q[0];
rz(3.13700389036024) q[0];
sx q[0];
rz(11.1176783800046) q[0];
rz(-0.121512420475483) q[2];
sx q[2];
rz(3.80282375414903) q[2];
sx q[2];
rz(11.4610550165097) q[2];
cx q[2],q[1];
rz(0.400075167417526) q[1];
sx q[1];
rz(4.14717522461946) q[1];
sx q[1];
rz(11.2838047504346) q[1];
rz(0.469623029232025) q[3];
sx q[3];
rz(4.22007385094697) q[3];
sx q[3];
rz(10.0520536661069) q[3];
cx q[3],q[2];
rz(1.23624706268311) q[2];
sx q[2];
rz(3.6422211249643) q[2];
sx q[2];
rz(9.53510366975471) q[2];
rz(0.825354397296906) q[3];
sx q[3];
rz(2.3767410834604) q[3];
sx q[3];
rz(10.1650300979535) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.16846191883087) q[0];
sx q[0];
rz(3.35947334964807) q[0];
sx q[0];
rz(10.3067807912748) q[0];
rz(-1.02856707572937) q[1];
sx q[1];
rz(2.587963732081) q[1];
sx q[1];
rz(11.6539547204892) q[1];
cx q[1],q[0];
rz(-0.00936326570808887) q[0];
sx q[0];
rz(3.5991266985708) q[0];
sx q[0];
rz(9.26079784928962) q[0];
rz(-0.120435111224651) q[2];
sx q[2];
rz(3.99748799403245) q[2];
sx q[2];
rz(9.46288409679338) q[2];
cx q[2],q[1];
rz(-0.470570236444473) q[1];
sx q[1];
rz(3.81524130900437) q[1];
sx q[1];
rz(8.87646911143466) q[1];
rz(0.498982042074203) q[3];
sx q[3];
rz(3.64078882535035) q[3];
sx q[3];
rz(10.04933152198) q[3];
cx q[3],q[2];
rz(2.26221108436584) q[2];
sx q[2];
rz(3.38962720532949) q[2];
sx q[2];
rz(8.60652718543216) q[2];
rz(0.360005587339401) q[3];
sx q[3];
rz(4.44029048283631) q[3];
sx q[3];
rz(9.46258757486149) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.705719470977783) q[0];
sx q[0];
rz(4.09214028914506) q[0];
sx q[0];
rz(10.5856202602307) q[0];
rz(0.970876514911652) q[1];
sx q[1];
rz(4.05738070805604) q[1];
sx q[1];
rz(9.53825574218437) q[1];
cx q[1],q[0];
rz(0.161008536815643) q[0];
sx q[0];
rz(4.11797395546968) q[0];
sx q[0];
rz(9.91765967606708) q[0];
rz(0.356747359037399) q[2];
sx q[2];
rz(3.10380195279653) q[2];
sx q[2];
rz(10.4228849172513) q[2];
cx q[2],q[1];
rz(1.74746882915497) q[1];
sx q[1];
rz(4.18167665799195) q[1];
sx q[1];
rz(10.2323540210645) q[1];
rz(-0.640796542167664) q[3];
sx q[3];
rz(4.17916813691194) q[3];
sx q[3];
rz(9.81742430328532) q[3];
cx q[3],q[2];
rz(0.111614413559437) q[2];
sx q[2];
rz(4.47776416142518) q[2];
sx q[2];
rz(10.3327769994657) q[2];
rz(0.326344311237335) q[3];
sx q[3];
rz(3.69644531806047) q[3];
sx q[3];
rz(9.84659100174114) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.570807635784149) q[0];
sx q[0];
rz(3.90358481009538) q[0];
sx q[0];
rz(10.4417814969937) q[0];
rz(0.486558884382248) q[1];
sx q[1];
rz(4.16688814957673) q[1];
sx q[1];
rz(9.51973254083797) q[1];
cx q[1],q[0];
rz(0.217903777956963) q[0];
sx q[0];
rz(2.99605401058728) q[0];
sx q[0];
rz(10.1654579996984) q[0];
rz(0.21910659968853) q[2];
sx q[2];
rz(3.89987513621385) q[2];
sx q[2];
rz(10.1667423009793) q[2];
cx q[2],q[1];
rz(0.508521139621735) q[1];
sx q[1];
rz(4.43852940400178) q[1];
sx q[1];
rz(10.6349411964337) q[1];
rz(1.26035916805267) q[3];
sx q[3];
rz(4.41239550908143) q[3];
sx q[3];
rz(9.83481568693324) q[3];
cx q[3],q[2];
rz(0.9359370470047) q[2];
sx q[2];
rz(3.7762452681833) q[2];
sx q[2];
rz(9.89457160829707) q[2];
rz(0.698224067687988) q[3];
sx q[3];
rz(5.05687383015687) q[3];
sx q[3];
rz(9.74495939015552) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.532233774662018) q[0];
sx q[0];
rz(2.40257522662217) q[0];
sx q[0];
rz(9.63960174321338) q[0];
rz(-0.23669445514679) q[1];
sx q[1];
rz(3.71859076817567) q[1];
sx q[1];
rz(9.844066476814) q[1];
cx q[1],q[0];
rz(0.400601655244827) q[0];
sx q[0];
rz(3.23232750047977) q[0];
sx q[0];
rz(8.89906094073459) q[0];
rz(1.88168787956238) q[2];
sx q[2];
rz(4.19294801552827) q[2];
sx q[2];
rz(8.63294098376437) q[2];
cx q[2],q[1];
rz(0.579455733299255) q[1];
sx q[1];
rz(3.53295028408105) q[1];
sx q[1];
rz(9.98387191294833) q[1];
rz(-0.198344424366951) q[3];
sx q[3];
rz(2.2738278230005) q[3];
sx q[3];
rz(11.2445994377057) q[3];
cx q[3],q[2];
rz(0.995701670646667) q[2];
sx q[2];
rz(3.27708430786664) q[2];
sx q[2];
rz(9.51616936772271) q[2];
rz(0.351783007383347) q[3];
sx q[3];
rz(3.88889035780961) q[3];
sx q[3];
rz(8.959564959995) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.92129373550415) q[0];
sx q[0];
rz(4.31381729443605) q[0];
sx q[0];
rz(10.6002909898679) q[0];
rz(0.579188644886017) q[1];
sx q[1];
rz(3.94828221400315) q[1];
sx q[1];
rz(9.61434141396686) q[1];
cx q[1],q[0];
rz(1.26111471652985) q[0];
sx q[0];
rz(4.91312554677064) q[0];
sx q[0];
rz(9.52331984638377) q[0];
rz(0.102645561099052) q[2];
sx q[2];
rz(2.43754390080506) q[2];
sx q[2];
rz(11.2789801120679) q[2];
cx q[2],q[1];
rz(0.71566641330719) q[1];
sx q[1];
rz(3.95968422492082) q[1];
sx q[1];
rz(10.9851500749509) q[1];
rz(0.104699693620205) q[3];
sx q[3];
rz(2.61879131396348) q[3];
sx q[3];
rz(10.3308498024861) q[3];
cx q[3],q[2];
rz(0.243894800543785) q[2];
sx q[2];
rz(4.44111851056153) q[2];
sx q[2];
rz(9.93068698643848) q[2];
rz(0.778877437114716) q[3];
sx q[3];
rz(3.53770309885079) q[3];
sx q[3];
rz(10.7194322109143) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.270957231521606) q[0];
sx q[0];
rz(4.19356182416017) q[0];
sx q[0];
rz(10.58027110099) q[0];
rz(-0.00118601543363184) q[1];
sx q[1];
rz(3.92918190558488) q[1];
sx q[1];
rz(9.82996413706943) q[1];
cx q[1],q[0];
rz(0.74086606502533) q[0];
sx q[0];
rz(2.40805456240708) q[0];
sx q[0];
rz(10.2597784757535) q[0];
rz(-0.343339323997498) q[2];
sx q[2];
rz(3.57309588988359) q[2];
sx q[2];
rz(10.4017020225446) q[2];
cx q[2],q[1];
rz(0.894541382789612) q[1];
sx q[1];
rz(3.90905848343904) q[1];
sx q[1];
rz(9.53453150986835) q[1];
rz(1.08069813251495) q[3];
sx q[3];
rz(2.99672720034654) q[3];
sx q[3];
rz(9.64986096917793) q[3];
cx q[3],q[2];
rz(0.817951679229736) q[2];
sx q[2];
rz(4.30887821515138) q[2];
sx q[2];
rz(9.94636384247943) q[2];
rz(-0.0966054350137711) q[3];
sx q[3];
rz(4.16028264363343) q[3];
sx q[3];
rz(9.91560933589145) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.38799649477005) q[0];
sx q[0];
rz(4.02194455464418) q[0];
sx q[0];
rz(9.53943592905208) q[0];
rz(1.82785868644714) q[1];
sx q[1];
rz(4.82721713383729) q[1];
sx q[1];
rz(9.55984455942317) q[1];
cx q[1],q[0];
rz(1.01975214481354) q[0];
sx q[0];
rz(4.49603560765321) q[0];
sx q[0];
rz(10.7309728622357) q[0];
rz(2.15536403656006) q[2];
sx q[2];
rz(3.60167142947251) q[2];
sx q[2];
rz(9.73420888780757) q[2];
cx q[2],q[1];
rz(1.9070930480957) q[1];
sx q[1];
rz(3.78428253729875) q[1];
sx q[1];
rz(9.69729698299571) q[1];
rz(-0.11718312650919) q[3];
sx q[3];
rz(2.57962671120698) q[3];
sx q[3];
rz(9.19807054697677) q[3];
cx q[3],q[2];
rz(-0.315852850675583) q[2];
sx q[2];
rz(4.35704723198945) q[2];
sx q[2];
rz(10.3703534960668) q[2];
rz(0.259975075721741) q[3];
sx q[3];
rz(4.16444042523439) q[3];
sx q[3];
rz(10.5528356790464) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.05396997928619) q[0];
sx q[0];
rz(4.23183432419831) q[0];
sx q[0];
rz(10.8277431487958) q[0];
rz(0.915790915489197) q[1];
sx q[1];
rz(3.75800153811509) q[1];
sx q[1];
rz(10.3403299212377) q[1];
cx q[1],q[0];
rz(-0.212857350707054) q[0];
sx q[0];
rz(2.97862507601316) q[0];
sx q[0];
rz(8.78446451424762) q[0];
rz(0.327387601137161) q[2];
sx q[2];
rz(3.6871579011255) q[2];
sx q[2];
rz(8.69395265578433) q[2];
cx q[2],q[1];
rz(1.82322227954865) q[1];
sx q[1];
rz(3.75311932166154) q[1];
sx q[1];
rz(8.67728844880267) q[1];
rz(1.40110182762146) q[3];
sx q[3];
rz(2.27329191763932) q[3];
sx q[3];
rz(9.19965738653346) q[3];
cx q[3],q[2];
rz(-0.123306840658188) q[2];
sx q[2];
rz(3.78458896477754) q[2];
sx q[2];
rz(9.0221590757291) q[2];
rz(-1.17665243148804) q[3];
sx q[3];
rz(3.42824965913827) q[3];
sx q[3];
rz(10.1927463769834) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.313189327716827) q[0];
sx q[0];
rz(3.93149110873277) q[0];
sx q[0];
rz(10.3405120730321) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(-0.340436518192291) q[1];
sx q[1];
rz(3.93991831143434) q[1];
sx q[1];
rz(10.7481257676999) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(0.2764752805233) q[2];
sx q[2];
rz(4.14294341404969) q[2];
sx q[2];
rz(10.2793879270475) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(0.798740327358246) q[3];
sx q[3];
rz(3.83161416848237) q[3];
sx q[3];
rz(10.5598190784375) q[3];
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
