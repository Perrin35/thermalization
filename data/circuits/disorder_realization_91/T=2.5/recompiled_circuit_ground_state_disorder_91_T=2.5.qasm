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
rz(-1.81023263931274) q[0];
sx q[0];
rz(4.39637640317018) q[0];
sx q[0];
rz(8.77883610724612) q[0];
rz(1.68260014057159) q[1];
sx q[1];
rz(3.26972200174863) q[1];
sx q[1];
rz(10.0929565787236) q[1];
cx q[1],q[0];
rz(0.0543391518294811) q[0];
sx q[0];
rz(2.87326604326303) q[0];
sx q[0];
rz(7.41560933589145) q[0];
rz(1.24348843097687) q[2];
sx q[2];
rz(7.23282256920869) q[2];
sx q[2];
rz(3.60840890406772) q[2];
cx q[2],q[1];
rz(-2.97578024864197) q[1];
sx q[1];
rz(6.96271673043305) q[1];
sx q[1];
rz(6.94289062022373) q[1];
rz(5.5623779296875) q[3];
sx q[3];
rz(1.94879833062226) q[3];
sx q[3];
rz(7.57021615504428) q[3];
cx q[3],q[2];
rz(-2.62319231033325) q[2];
sx q[2];
rz(3.79700848658616) q[2];
sx q[2];
rz(8.34168121813937) q[2];
rz(-3.3969349861145) q[3];
sx q[3];
rz(5.3861929496103) q[3];
sx q[3];
rz(13.7258796453397) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.97936797142029) q[0];
sx q[0];
rz(3.00016820629174) q[0];
sx q[0];
rz(9.39119957610174) q[0];
rz(1.13831877708435) q[1];
sx q[1];
rz(4.13616290886933) q[1];
sx q[1];
rz(12.8033270597379) q[1];
cx q[1],q[0];
rz(0.367192715406418) q[0];
sx q[0];
rz(6.01472345192964) q[0];
sx q[0];
rz(11.5422572851102) q[0];
rz(0.476322501897812) q[2];
sx q[2];
rz(4.54541316826875) q[2];
sx q[2];
rz(10.1379267930905) q[2];
cx q[2],q[1];
rz(-1.10190558433533) q[1];
sx q[1];
rz(6.47727742989595) q[1];
sx q[1];
rz(9.55691685377761) q[1];
rz(3.77443885803223) q[3];
sx q[3];
rz(2.45996353228623) q[3];
sx q[3];
rz(11.0107811450879) q[3];
cx q[3],q[2];
rz(2.52915501594543) q[2];
sx q[2];
rz(2.09450224240357) q[2];
sx q[2];
rz(9.30442656426832) q[2];
rz(-0.585931658744812) q[3];
sx q[3];
rz(4.52208581765229) q[3];
sx q[3];
rz(11.2980414390485) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.8782069683075) q[0];
sx q[0];
rz(4.10739121039445) q[0];
sx q[0];
rz(10.4129218220632) q[0];
rz(1.49096095561981) q[1];
sx q[1];
rz(5.56628027756745) q[1];
sx q[1];
rz(7.83685538767978) q[1];
cx q[1],q[0];
rz(-0.40000382065773) q[0];
sx q[0];
rz(3.15753879223997) q[0];
sx q[0];
rz(11.0647072553556) q[0];
rz(-1.70519804954529) q[2];
sx q[2];
rz(4.64407947857911) q[2];
sx q[2];
rz(10.8022843360822) q[2];
cx q[2],q[1];
rz(-2.88042593002319) q[1];
sx q[1];
rz(6.0203813632303) q[1];
sx q[1];
rz(9.31402397751018) q[1];
rz(-0.421093493700027) q[3];
sx q[3];
rz(5.23947969277436) q[3];
sx q[3];
rz(8.52839622496768) q[3];
cx q[3],q[2];
rz(-8.58738708496094) q[2];
sx q[2];
rz(1.34342554410035) q[2];
sx q[2];
rz(12.4957911729734) q[2];
rz(-0.48921674489975) q[3];
sx q[3];
rz(4.13240465720231) q[3];
sx q[3];
rz(12.418925023071) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.95511424541473) q[0];
sx q[0];
rz(5.58772793610627) q[0];
sx q[0];
rz(7.68408463000461) q[0];
rz(3.31162285804749) q[1];
sx q[1];
rz(1.73151198227937) q[1];
sx q[1];
rz(8.19164285659) q[1];
cx q[1],q[0];
rz(-0.296516686677933) q[0];
sx q[0];
rz(4.29238465626771) q[0];
sx q[0];
rz(10.8588909864347) q[0];
rz(-0.581572115421295) q[2];
sx q[2];
rz(4.08245119650895) q[2];
sx q[2];
rz(9.41897083808809) q[2];
cx q[2],q[1];
rz(8.48272132873535) q[1];
sx q[1];
rz(2.8615073283487) q[1];
sx q[1];
rz(6.33732078074619) q[1];
rz(-1.1063152551651) q[3];
sx q[3];
rz(4.72314503987367) q[3];
sx q[3];
rz(10.6841983556668) q[3];
cx q[3],q[2];
rz(7.082923412323) q[2];
sx q[2];
rz(1.43504718144471) q[2];
sx q[2];
rz(5.65356872080966) q[2];
rz(-0.136842116713524) q[3];
sx q[3];
rz(5.65356579621369) q[3];
sx q[3];
rz(10.8401661872785) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.283023089170456) q[0];
sx q[0];
rz(3.6498254259401) q[0];
sx q[0];
rz(9.3224757745783) q[0];
rz(1.49818587303162) q[1];
sx q[1];
rz(1.30032769043977) q[1];
sx q[1];
rz(12.209265446655) q[1];
cx q[1],q[0];
rz(1.50955605506897) q[0];
sx q[0];
rz(0.311838301020213) q[0];
sx q[0];
rz(10.2300708651464) q[0];
rz(1.35499143600464) q[2];
sx q[2];
rz(3.26769548852975) q[2];
sx q[2];
rz(9.4377735465686) q[2];
cx q[2],q[1];
rz(-2.63729238510132) q[1];
sx q[1];
rz(5.01070395310456) q[1];
sx q[1];
rz(8.38506577014133) q[1];
rz(2.65642046928406) q[3];
sx q[3];
rz(1.96868041356141) q[3];
sx q[3];
rz(11.4317717313687) q[3];
cx q[3],q[2];
rz(2.66852784156799) q[2];
sx q[2];
rz(1.7230955680185) q[2];
sx q[2];
rz(14.3853745222013) q[2];
rz(1.43329167366028) q[3];
sx q[3];
rz(4.96634105046327) q[3];
sx q[3];
rz(9.58874455689594) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.543561637401581) q[0];
sx q[0];
rz(4.14154431422288) q[0];
sx q[0];
rz(7.7213823556821) q[0];
rz(4.3819055557251) q[1];
sx q[1];
rz(3.79642436106736) q[1];
sx q[1];
rz(13.9192194700162) q[1];
cx q[1],q[0];
rz(0.082093246281147) q[0];
sx q[0];
rz(1.66498604615266) q[0];
sx q[0];
rz(12.159417128555) q[0];
rz(-1.04463052749634) q[2];
sx q[2];
rz(4.858502419787) q[2];
sx q[2];
rz(9.86175472139522) q[2];
cx q[2],q[1];
rz(-1.28619146347046) q[1];
sx q[1];
rz(2.07780650456483) q[1];
sx q[1];
rz(11.4195095062177) q[1];
rz(1.30451726913452) q[3];
sx q[3];
rz(3.54972800810868) q[3];
sx q[3];
rz(8.54165402650043) q[3];
cx q[3],q[2];
rz(-1.21443402767181) q[2];
sx q[2];
rz(2.77756014664704) q[2];
sx q[2];
rz(6.67361710070773) q[2];
rz(-1.82921421527863) q[3];
sx q[3];
rz(3.95268419583375) q[3];
sx q[3];
rz(10.2446027159612) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.445359617471695) q[0];
sx q[0];
rz(4.06548252900178) q[0];
sx q[0];
rz(7.69611415862247) q[0];
rz(-1.67848086357117) q[1];
sx q[1];
rz(4.69722727139527) q[1];
sx q[1];
rz(11.9268371820371) q[1];
cx q[1],q[0];
rz(0.937767505645752) q[0];
sx q[0];
rz(7.26221719582612) q[0];
sx q[0];
rz(10.5401381015699) q[0];
rz(3.06601428985596) q[2];
sx q[2];
rz(4.86267081101472) q[2];
sx q[2];
rz(7.7102328300397) q[2];
cx q[2],q[1];
rz(-2.8486316204071) q[1];
sx q[1];
rz(5.34103408654267) q[1];
sx q[1];
rz(12.8958666086118) q[1];
rz(0.280847132205963) q[3];
sx q[3];
rz(5.87434712250764) q[3];
sx q[3];
rz(10.5146896600644) q[3];
cx q[3],q[2];
rz(2.92452359199524) q[2];
sx q[2];
rz(4.4704419692331) q[2];
sx q[2];
rz(9.8769763469617) q[2];
rz(0.222981199622154) q[3];
sx q[3];
rz(3.42295134265954) q[3];
sx q[3];
rz(12.7217173337857) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(3.01057910919189) q[0];
sx q[0];
rz(5.36088720162446) q[0];
sx q[0];
rz(9.58670122026607) q[0];
rz(2.79363918304443) q[1];
sx q[1];
rz(4.41569248040254) q[1];
sx q[1];
rz(15.4735698461454) q[1];
cx q[1],q[0];
rz(2.16805553436279) q[0];
sx q[0];
rz(4.13482031424577) q[0];
sx q[0];
rz(10.319080746166) q[0];
rz(-1.21974813938141) q[2];
sx q[2];
rz(3.89294150670106) q[2];
sx q[2];
rz(10.0883273839872) q[2];
cx q[2],q[1];
rz(1.76633155345917) q[1];
sx q[1];
rz(5.34577074845368) q[1];
sx q[1];
rz(13.9042768239896) q[1];
rz(1.33998811244965) q[3];
sx q[3];
rz(2.85770204861695) q[3];
sx q[3];
rz(9.89764512180492) q[3];
cx q[3],q[2];
rz(-0.286354124546051) q[2];
sx q[2];
rz(5.21205607255036) q[2];
sx q[2];
rz(11.9422428369443) q[2];
rz(0.74328076839447) q[3];
sx q[3];
rz(5.28947654564912) q[3];
sx q[3];
rz(10.1144189000051) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.2766991853714) q[0];
sx q[0];
rz(4.58692160447175) q[0];
sx q[0];
rz(7.37010190486118) q[0];
rz(-1.90896058082581) q[1];
sx q[1];
rz(7.39092126687104) q[1];
sx q[1];
rz(4.98082015513583) q[1];
cx q[1],q[0];
rz(0.248274579644203) q[0];
sx q[0];
rz(3.38015462656552) q[0];
sx q[0];
rz(11.493625855438) q[0];
rz(1.21643054485321) q[2];
sx q[2];
rz(5.27949634392793) q[2];
sx q[2];
rz(10.7504156589429) q[2];
cx q[2],q[1];
rz(1.2375396490097) q[1];
sx q[1];
rz(4.31892052491242) q[1];
sx q[1];
rz(11.572254395477) q[1];
rz(1.59713852405548) q[3];
sx q[3];
rz(3.76479366620118) q[3];
sx q[3];
rz(10.6818105935971) q[3];
cx q[3],q[2];
rz(2.67526006698608) q[2];
sx q[2];
rz(4.39753547509248) q[2];
sx q[2];
rz(9.01042509674236) q[2];
rz(0.770533442497253) q[3];
sx q[3];
rz(4.56376782258088) q[3];
sx q[3];
rz(8.49111548661395) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.829380214214325) q[0];
sx q[0];
rz(0.820136221247264) q[0];
sx q[0];
rz(9.89341548680469) q[0];
rz(4.41523694992065) q[1];
sx q[1];
rz(8.14866271813447) q[1];
sx q[1];
rz(12.2548849344175) q[1];
cx q[1],q[0];
rz(3.71989345550537) q[0];
sx q[0];
rz(1.51307633717591) q[0];
sx q[0];
rz(7.00194237231418) q[0];
rz(0.0869290977716446) q[2];
sx q[2];
rz(4.82918253739411) q[2];
sx q[2];
rz(11.6051139593045) q[2];
cx q[2],q[1];
rz(-1.12324166297913) q[1];
sx q[1];
rz(9.91453281243379) q[1];
sx q[1];
rz(9.32547452150985) q[1];
rz(-0.488481611013412) q[3];
sx q[3];
rz(2.04601076443727) q[3];
sx q[3];
rz(12.1679653882901) q[3];
cx q[3],q[2];
rz(1.10625386238098) q[2];
sx q[2];
rz(4.33572748501832) q[2];
sx q[2];
rz(11.1733250379483) q[2];
rz(2.2140908241272) q[3];
sx q[3];
rz(4.71908834775025) q[3];
sx q[3];
rz(10.2849910020749) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.37921118736267) q[0];
sx q[0];
rz(5.03978279431398) q[0];
sx q[0];
rz(10.9671963214795) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(-2.42770910263062) q[1];
sx q[1];
rz(-0.857029525441579) q[1];
sx q[1];
rz(11.6804175138394) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(2.91412901878357) q[2];
sx q[2];
rz(4.59039655526216) q[2];
sx q[2];
rz(9.95753655432864) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(-1.0925532579422) q[3];
sx q[3];
rz(0.577742012339183) q[3];
sx q[3];
rz(7.74674043654605) q[3];
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
