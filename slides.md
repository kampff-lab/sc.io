---
layout: slides
title: Slides
permalink: /slides/
---

{% for slide in site.intro %}
<section>
<h3>{{ slide.title }}</h3>
{% if slide.image %}
<div style="width:100%">
  <div style="width:60%;float:left">{{ slide.content }}</div>
  <div style="width:40%;float:right"><img src="{{ slide.image }}" /></div>
</div>
{% else %}
{{ slide.content }}
{% endif %}
</section>
{% endfor %}
